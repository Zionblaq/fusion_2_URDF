"""
Fusion_to_URDF
==============
Exports the active Fusion design to a URDF package (urdf + STL meshes).

Flow:
  1. CHECK  - reads the whole design and lists everything that would break a URDF,
              each with the exact fix. If there are errors, nothing is exported.
  2. EXPORT - if the design passes, it writes the URDF and meshes.
  3. VERIFY - re-reads the written URDF and checks it.
  4. If anything fails at any stage, a report says which stage / which part failed,
     why that usually happens, and what to correct.

It NEVER modifies your design (no renaming, no copying components).
Units: Fusion API = cm / kg. URDF = metres / kg.
"""
import adsk.core, adsk.fusion, traceback, re, os, math, struct, datetime
from collections import deque
from xml.etree import ElementTree as ET
from xml.dom import minidom

# ----------------------------- settings -----------------------------
MESH_QUALITY = 'Low'        # 'Low', 'Medium' or 'High'
DEFAULT_EFFORT = 100.0
DEFAULT_VELOCITY = 100.0
REPORT_PATH = os.path.join(os.path.expanduser('~'), 'Documents', 'urdf_errors.txt')

JOINT_NAMES = {0: 'Rigid', 1: 'Revolute', 2: 'Slider', 3: 'Cylindrical',
               4: 'Pin-Slot', 5: 'Planar', 6: 'Ball'}
SUPPORTED = (0, 1, 2)


# ----------------------------- helpers ------------------------------
class Report:
    def __init__(self):
        self.errors, self.warnings, self.infos = [], [], []
        self.tb = ''

    def error(self, problem, fix):
        self.errors.append((problem, fix))

    def warn(self, problem, fix=''):
        self.warnings.append((problem, fix))

    def info(self, text):
        self.infos.append(text)

    def text(self, doc, headline):
        L = ['=' * 66,
             'FUSION -> URDF   %s' % datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
             'Document: %s' % doc, '=' * 66]
        if self.errors:
            L.append('\nALL %d ERRORS AT A GLANCE (fix them all, then run again):' % len(self.errors))
            for i, (p, f) in enumerate(self.errors, 1):
                L.append('  [ ] %d. %s' % (i, p if len(p) < 110 else p[:107] + '...'))
            L.append('\nERRORS IN DETAIL: %d' % len(self.errors))
            for i, (p, f) in enumerate(self.errors, 1):
                L.append('\n%d. %s\n   -> DO THIS: %s' % (i, p, f))
        if self.warnings:
            L.append('\nWARNINGS (export continues, but check): %d' % len(self.warnings))
            for i, (p, f) in enumerate(self.warnings, 1):
                L.append('\n%d. %s' % (i, p) + ('\n   -> %s' % f if f else ''))
        if self.infos:
            L.append('\nNOTES:')
            for t in self.infos:
                L.append(' - ' + t)
        if self.tb:
            L.append('\nTECHNICAL DETAIL (send this if you need help):\n' + self.tb)
        L += ['\n' + '-' * 66, headline, '-' * 66]
        return '\n'.join(L)


def sanitize(s):
    return re.sub(r'[^A-Za-z0-9_]+', '_', s).strip('_') or 'unnamed'


def fmt(v):
    v = float(v)
    if abs(v) < 1e-12:
        v = 0.0
    return '%.8g' % v


def vec(v):
    return ' '.join(fmt(x) for x in v)


def has_geometry(occ):
    """True if the occurrence has a body itself or anywhere in its nested parts."""
    try:
        if occ.bRepBodies.count > 0:
            return True
        for c in occ.childOccurrences:
            if has_geometry(c):
                return True
    except:
        pass
    return False


def top_of(occ):
    while occ.assemblyContext is not None:
        occ = occ.assemblyContext
    return occ


def joint_origin_cm(j):
    g = j.geometryOrOriginOne
    jo = adsk.fusion.JointOrigin.cast(g)
    if jo:
        return jo.geometry.origin
    jg = adsk.fusion.JointGeometry.cast(g)
    if jg:
        return jg.origin
    raise ValueError('unknown joint geometry type')


def joint_axis(j, jt):
    jm = j.jointMotion
    if jt == 1:
        v = adsk.fusion.RevoluteJointMotion.cast(jm).rotationAxisVector
    else:
        v = adsk.fusion.SliderJointMotion.cast(jm).slideDirectionVector
    n = math.sqrt(v.x ** 2 + v.y ** 2 + v.z ** 2) or 1.0
    return (v.x / n, v.y / n, v.z / n)


# ----------------------------- STL helpers --------------------------
def stl_bounds(path):
    with open(path, 'rb') as f:
        data = f.read()
    mn = [1e30] * 3
    mx = [-1e30] * 3
    if len(data) >= 84:
        n = struct.unpack_from('<I', data, 80)[0]
        if len(data) == 84 + 50 * n and n > 0:
            for rec in struct.iter_unpack('<12fH', data[84:84 + 50 * n]):
                for k in range(3):
                    a, b, c = rec[3 + k], rec[6 + k], rec[9 + k]
                    lo, hi = min(a, b, c), max(a, b, c)
                    if lo < mn[k]:
                        mn[k] = lo
                    if hi > mx[k]:
                        mx[k] = hi
            return mn, mx
    # ASCII fallback
    txt = data.decode('latin-1', 'ignore')
    for m in re.finditer(r'vertex\s+(\S+)\s+(\S+)\s+(\S+)', txt):
        for k in range(3):
            x = float(m.group(k + 1))
            mn[k] = min(mn[k], x)
            mx[k] = max(mx[k], x)
    if mn[0] > mx[0]:
        raise ValueError('STL file is empty')
    return mn, mx


def rpy_from_matrix(R):
    pitch = -math.asin(max(-1.0, min(1.0, R[2][0])))
    if abs(math.cos(pitch)) > 1e-6:
        roll = math.atan2(R[2][1], R[2][2])
        yaw = math.atan2(R[1][0], R[0][0])
    else:
        roll = 0.0
        yaw = math.atan2(-R[0][1], R[1][1])
    return (roll, pitch, yaw)


def mesh_placement(stl_path, occ, link_world):
    """
    Returns (scale_to_metres, visual_xyz, visual_rpy, mode)
    mode: 'world'    mesh is already in assembly coordinates
          'local'    mesh is in the part's own frame; placed with the part's transform (position + rotation)
          'fallback' could not be matched; re-centred by bounding box (rotation not handled)
    """
    bb = occ.boundingBox
    wmin = [bb.minPoint.x * 0.01, bb.minPoint.y * 0.01, bb.minPoint.z * 0.01]
    wmax = [bb.maxPoint.x * 0.01, bb.maxPoint.y * 0.01, bb.maxPoint.z * 0.01]
    ext_w = max(wmax[i] - wmin[i] for i in range(3))
    smin, smax = stl_bounds(stl_path)
    ext_s = max(smax[i] - smin[i] for i in range(3))
    if ext_s <= 0:
        raise ValueError('mesh has zero size')
    if ext_w > 0:
        ratio = ext_w / ext_s
        scale = min([1.0, 0.1, 0.01, 0.001], key=lambda s: abs(math.log(ratio / s)))
    else:
        scale = 0.001
    tol = 0.03 * ext_w + 0.002

    def close(lo, hi):
        return all(abs(lo[i] - wmin[i]) <= tol and abs(hi[i] - wmax[i]) <= tol for i in range(3))

    # 1) mesh already in world coordinates?
    lo = [smin[i] * scale for i in range(3)]
    hi = [smax[i] * scale for i in range(3)]
    if close(lo, hi):
        return scale, [-link_world[i] for i in range(3)], [0.0, 0.0, 0.0], 'world'

    # 2) mesh in the part's own frame: world = R * local + t
    T = occ.transform
    R = [[T.getCell(r, c) for c in range(3)] for r in range(3)]
    t = [T.getCell(r, 3) * 0.01 for r in range(3)]
    xs = (smin[0] * scale, smax[0] * scale)
    ys = (smin[1] * scale, smax[1] * scale)
    zs = (smin[2] * scale, smax[2] * scale)
    pts = [[R[r][0] * x + R[r][1] * y + R[r][2] * z + t[r] for r in range(3)]
           for x in xs for y in ys for z in zs]
    lo2 = [min(p[i] for p in pts) for i in range(3)]
    hi2 = [max(p[i] for p in pts) for i in range(3)]
    if close(lo2, hi2):
        return scale, [t[i] - link_world[i] for i in range(3)], list(rpy_from_matrix(R)), 'local'

    # 3) fallback: centre by bounding box (position only)
    mc = [(smin[i] + smax[i]) / 2.0 * scale for i in range(3)]
    wc = [(wmin[i] + wmax[i]) / 2.0 for i in range(3)]
    return scale, [wc[i] - link_world[i] - mc[i] for i in range(3)], [0.0, 0.0, 0.0], 'fallback'


# ----------------------------- STAGE 1: CHECK ------------------------
def analyze(design, rep, st):
    root = design.rootComponent
    st['stage'], st['item'] = 'checking the design', ''
    top = list(root.occurrences)
    cand = [o for o in top if has_geometry(o)]

    if not cand:
        rep.error('The design has no components that contain bodies.',
                  'Create components with solid bodies at the top level of the Browser.')
        return None

    # leftovers from the old exporter
    for o in top:
        cn = o.component.name
        if cn == 'old_component' or o.name.lower().startswith('base_link (') or \
                (re.match(r'^Component\d+$', cn) and not has_geometry(o)):
            rep.warn('"%s" looks like a leftover from a failed export with the old script.' % o.name,
                     'Undo (Ctrl+Z) or reopen your saved copy, so you do not export junk.')

    # base link
    bl = [o for o in cand if o.component.name == 'base_link']
    base = None
    if len(bl) > 1:
        rep.error('More than one component is named base_link: %s.' % [o.name for o in bl],
                  'Keep one component named base_link and rename the others.')
    elif len(bl) == 1:
        base = bl[0]
    else:
        g = [o for o in cand if o.isGrounded]
        if len(g) == 1:
            base = g[0]
            rep.warn('No component is named base_link. Using the grounded component "%s" as the base.' % base.name,
                     'To silence this, rename that component to base_link.')
        else:
            rep.error('No component is named base_link, and I cannot tell which one is the base '
                      '(%d grounded components).' % len(g),
                      'Rename your main fixed component (the chassis) to base_link, and Ground it '
                      '(right-click > Ground).')
    # link names (continue checking even if base_link is missing)
    names = {}
    used = set()
    if base is not None:
        names[base.name] = 'base_link'
        used.add('base_link')
    for o in cand:
        if o.name in names:
            continue
        n = sanitize(o.name)
        b, k = n, 2
        while n in used:
            n = '%s_%d' % (b, k)
            k += 1
        used.add(n)
        names[o.name] = n
        if n != o.name.replace(':', '_').replace(' ', '_'):
            pass
    renamed = [(o.name, names[o.name]) for o in cand if names[o.name] != o.name]
    if renamed:
        rep.info('Names cleaned for URDF: ' + ', '.join('%s -> %s' % r for r in renamed[:12]) +
                 (' ...' if len(renamed) > 12 else ''))

    # rigid groups / as-built / sub-component joints
    try:
        if root.rigidGroups.count > 0:
            rep.error('%d Rigid Group(s) found. They are not joints, so the parts end up disconnected.'
                      % root.rigidGroups.count,
                      'Delete each Rigid Group and add a Rigid joint between the two components.')
    except:
        pass
    try:
        if root.asBuiltJoints.count > 0:
            rep.warn('%d As-built joint(s) found. They are ignored by this exporter.'
                     % root.asBuiltJoints.count,
                     'Recreate them with Assemble > Joint if the parts must be connected.')
    except:
        pass
    try:
        sub = sum(c.joints.count for c in design.allComponents if c != root)
        if sub:
            rep.warn('%d joint(s) live inside sub-components and are ignored.' % sub,
                     'Recreate them in the top-level design.')
    except:
        pass

    # joints
    joints = []
    unresolved = []
    seen_pairs = {}
    st['stage'] = 'reading joints'
    for j in root.joints:
        st['item'] = j.name
        try:
            a, b = j.occurrenceOne, j.occurrenceTwo
            jt = int(j.jointMotion.jointType)
        except:
            rep.error('Joint "%s" cannot be read.' % j.name, 'Delete it and recreate it.')
            unresolved.append(j.name)
            continue
        if a is None or b is None:
            rep.error('Joint "%s" connects to the origin instead of a component.' % j.name,
                      'Delete it and recreate it between two components.')
            unresolved.append(j.name)
            continue
        ta, tb = top_of(a), top_of(b)
        bad = [t for t in (ta, tb) if t.name not in names]
        if bad:
            rep.error('Joint "%s" uses "%s", which has no bodies anywhere inside it.' % (j.name, bad[0].name),
                      'Put a body in that component, or joint a different component.')
            unresolved.append(j.name)
            continue
        if a.assemblyContext is not None or b.assemblyContext is not None:
            rep.warn('Joint "%s" is attached to a nested part; I treated it as a joint on the '
                     'top-level component (%s <-> %s).' % (j.name, ta.name, tb.name))
        la, lb = names[ta.name], names[tb.name]
        if la == lb:
            rep.warn('Joint "%s" connects a component to itself and was ignored.' % j.name)
            continue
        if jt not in SUPPORTED:
            rep.error('Joint "%s" is a %s joint, which URDF export here does not support.'
                      % (j.name, JOINT_NAMES.get(jt, jt)),
                      'Replace it with a Rigid, Revolute or Slider joint.')
        key = frozenset((la, lb))
        if key in seen_pairs:
            rep.warn('Joint "%s" duplicates "%s" (same two parts) and was ignored.'
                     % (j.name, seen_pairs[key]),
                     'Delete the extra joint in Fusion to keep things clean.')
            continue
        seen_pairs[key] = j.name

        lo = hi = None
        try:
            joint_origin_cm(j)
            if jt in (1, 2):
                joint_axis(j, jt)
        except Exception as ex:
            rep.error('Joint "%s": Fusion cannot give its position/axis (%s).' % (j.name, ex),
                      'Edit the joint and re-pick its geometry, or delete and recreate it.')
        try:
            if jt == 1:
                rl = adsk.fusion.RevoluteJointMotion.cast(j.jointMotion).rotationLimits
                if rl.isMinimumValueEnabled and rl.isMaximumValueEnabled:
                    lo, hi = rl.minimumValue, rl.maximumValue
                else:
                    rep.info('Revolute "%s" has no limits, exported as continuous '
                             '(fine for wheels; add Motion > Limits for hinges).' % j.name)
            elif jt == 2:
                sl = adsk.fusion.SliderJointMotion.cast(j.jointMotion).slideLimits
                if sl.isMinimumValueEnabled and sl.isMaximumValueEnabled:
                    lo, hi = sl.minimumValue * 0.01, sl.maximumValue * 0.01
                else:
                    rep.error('Slider joint "%s" has no limits. URDF prismatic joints need limits.' % j.name,
                              'Edit the joint, enable Motion > Limits and set minimum and maximum.')
        except:
            rep.warn('Could not read the limits of joint "%s".' % j.name)
        joints.append({'j': j, 'a': la, 'b': lb, 'jt': jt, 'lo': lo, 'hi': hi, 'name': j.name})

    # graph / tree
    note = ''
    if unresolved:
        note = (' (Some joints could not be used because of errors above, so this may disappear '
                'once those are fixed: %s)' % ', '.join(unresolved[:8]))
    adj = {n: [] for n in names.values()}
    for i, e in enumerate(joints):
        adj[e['a']].append((e['b'], i))
        adj[e['b']].append((e['a'], i))
    seen, order, tree, excluded = set(), [], [], []
    if base is not None:
        seen = {'base_link'}
        order = ['base_link']
        used_edges = set()
        dq = deque(['base_link'])
        while dq:
            u = dq.popleft()
            for v, i in adj[u]:
                if v not in seen:
                    seen.add(v)
                    order.append(v)
                    used_edges.add(i)
                    tree.append((u, v, i))
                    dq.append(v)
        for i, e in enumerate(joints):
            if i not in used_edges and e['a'] in seen and e['b'] in seen:
                rep.error('Joint "%s" closes a loop (%s <-> %s). URDF must be a tree, no loops.'
                          % (e['name'], e['a'], e['b']),
                          'Delete this joint (or break the loop elsewhere, e.g. 4-bar linkages).')
        for n, lst in adj.items():
            if n in seen:
                continue
            if lst:
                rep.error('"%s" has joints but is not connected to base_link.%s' % (n, note),
                          'Add a Rigid joint from "%s" (select it first) to a part that is connected '
                          'to base_link (select that second).' % n)
            else:
                excluded.append(n)
        if excluded:
            rep.warn('%d component(s) have no joints and were left OUT of the URDF: %s'
                     % (len(excluded), ', '.join(excluded)),
                     'To include one, add a Rigid joint from it to the part it is mounted on.')
        if not tree and not rep.errors:
            rep.error('No joint connects anything to base_link.',
                      'Add joints so every link is connected back to base_link.')
    else:
        groups, done = [], set()
        for n in adj:
            if n in done or not adj[n]:
                continue
            g, dq = {n}, deque([n])
            while dq:
                u = dq.popleft()
                for v, i in adj[u]:
                    if v not in g:
                        g.add(v)
                        dq.append(v)
            done |= g
            groups.append(sorted(g))
        rep.warn('Connectivity could not be fully checked because there is no base_link. '
                 'Connected groups right now: %s' % ('; '.join('[%s]' % ', '.join(g) for g in groups) or 'none'),
                 'Fix the base_link error first, then run again to check that every group connects to it.')
        return None

    rn = re.sub(r'\s+v\d+$', '', root.name.strip())
    robot_name = sanitize(rn) or 'robot'
    return {'root': root, 'occ': {names[o.name]: o for o in cand}, 'order': order,
            'tree': tree, 'joints': joints, 'robot': robot_name}


# ----------------------------- STAGE 2: EXPORT -----------------------
def export(design, model, outdir, rep, st):
    pkg = model['robot'].lower() + '_description'
    pdir = os.path.join(outdir, pkg)
    urdf_dir = os.path.join(pdir, 'urdf')
    mesh_dir = os.path.join(pdir, 'meshes')
    st['stage'], st['item'] = 'creating folders', pdir
    os.makedirs(urdf_dir, exist_ok=True)
    os.makedirs(mesh_dir, exist_ok=True)

    # link frames: base at world origin, others at their joint origin
    link_world = {'base_link': (0.0, 0.0, 0.0)}
    jinfo = {}
    st['stage'] = 'reading joint positions and axes'
    for parent, child, i in model['tree']:
        e = model['joints'][i]
        st['item'] = e['name']
        o = joint_origin_cm(e['j'])
        link_world[child] = (o.x * 0.01, o.y * 0.01, o.z * 0.01)
        flipped = (parent == e['a'])  # exporter convention: parent = occurrenceTwo
        axis = joint_axis(e['j'], e['jt']) if e['jt'] in (1, 2) else None
        lo, hi = e['lo'], e['hi']
        if flipped:
            if axis:
                axis = tuple(-x for x in axis)
            if lo is not None:
                lo, hi = -hi, -lo
        jinfo[child] = {'parent': parent, 'e': e, 'axis': axis, 'lo': lo, 'hi': hi}

    data = {}
    modes = {}
    em = design.exportManager
    q = {'Low': adsk.fusion.MeshRefinementSettings.MeshRefinementLow,
         'Medium': adsk.fusion.MeshRefinementSettings.MeshRefinementMedium,
         'High': adsk.fusion.MeshRefinementSettings.MeshRefinementHigh}[MESH_QUALITY]

    for name in model['order']:
        occ = model['occ'][name]
        lw = link_world[name]
        d = {}

        st['stage'], st['item'] = 'calculating mass and inertia', name
        try:
            pp = occ.getPhysicalProperties(adsk.fusion.CalculationAccuracy.MediumCalculationAccuracy)
            m = pp.mass
            c = pp.centerOfMass
            ok, xx, yy, zz, xy, yz, xz = pp.getXYZMomentsOfInertia()
            x, y, z = c.x, c.y, c.z
            ixx = (xx - m * (y * y + z * z)) * 1e-4
            iyy = (yy - m * (x * x + z * z)) * 1e-4
            izz = (zz - m * (x * x + y * y)) * 1e-4
            ixy = (xy + m * x * y) * 1e-4
            iyz = (yz + m * y * z) * 1e-4
            ixz = (xz + m * x * z) * 1e-4
            com = (x * 0.01, y * 0.01, z * 0.01)
            if m <= 0:
                raise ValueError('mass is zero')
            if ixx <= 0 or iyy <= 0 or izz <= 0:
                raise ValueError('inertia is not valid')
        except Exception as ex:
            m = 0.1
            com = (lw[0], lw[1], lw[2])
            ixx = iyy = izz = 1e-4
            ixy = iyz = ixz = 0.0
            rep.warn('Could not get a valid mass/inertia for "%s" (%s). Used 0.1 kg and a small '
                     'default inertia.' % (name, ex),
                     'Check the component has a solid body and a material. Right-click > Properties '
                     'should show a mass.')
        d.update(mass=m, com_rel=[com[i] - lw[i] for i in range(3)],
                 ixx=ixx, iyy=iyy, izz=izz, ixy=ixy, iyz=iyz, ixz=ixz)

        st['stage'], st['item'] = 'exporting the STL mesh', name
        try:
            path = os.path.join(mesh_dir, name + '.stl')
            if os.path.exists(path):
                os.remove(path)
            opts = em.createSTLExportOptions(occ, path)
            opts.sendToPrintUtility = False
            opts.isBinaryFormat = True
            opts.meshRefinement = q
            if not em.execute(opts):
                raise RuntimeError('Fusion refused to export the mesh')
            if not os.path.exists(path):
                alt = path[:-4]
                if os.path.exists(alt):
                    os.replace(alt, path)
                else:
                    raise RuntimeError('mesh file was not created')
            scale, vis, vrpy, mode = mesh_placement(path, occ, lw)
            d.update(scale=scale, vis=vis, vis_rpy=vrpy)
            modes.setdefault(mode, []).append(name)
            data[name] = d
        except Exception as ex:
            rep.error('The mesh for "%s" could not be exported (%s).' % (name, ex),
                      'Check "%s" has visible solid bodies. Try right-click the component > Save as '
                      'Mesh. If that fails too, repair the body (Modify > Repair Body) or simplify it.' % name)

    if rep.errors:
        return None
    if modes.get('local'):
        rep.info('%d mesh(es) were exported in each part\'s own frame; they were placed using the '
                 'part\'s position and rotation in the assembly.' % len(modes['local']))
    if modes.get('fallback'):
        rep.warn('The position of %d mesh(es) could not be matched to the assembly, so they were '
                 'centred by bounding box only (rotation ignored): %s'
                 % (len(modes['fallback']), ', '.join(modes['fallback'])),
                 'Open the URDF in a viewer. If these parts look rotated or misplaced, tell me and '
                 'send the report.')

    st['stage'], st['item'] = 'writing the URDF files', ''
    for suffix, prefix in (('', 'package://%s/meshes/' % pkg), ('_relative', '../meshes/')):
        xml = build_urdf(model, data, link_world, jinfo, prefix)
        with open(os.path.join(urdf_dir, model['robot'] + suffix + '.urdf'), 'w', encoding='utf-8') as f:
            f.write(xml)
    write_package_files(pdir, pkg)
    return pdir, os.path.join(urdf_dir, model['robot'] + '.urdf'), data, link_world


def build_urdf(model, data, link_world, jinfo, prefix):
    robot = ET.Element('robot', {'name': model['robot']})
    mat = ET.SubElement(robot, 'material', {'name': 'default_gray'})
    ET.SubElement(mat, 'color', {'rgba': '0.7 0.7 0.7 1.0'})
    for name in model['order']:
        d = data[name]
        link = ET.SubElement(robot, 'link', {'name': name})
        inr = ET.SubElement(link, 'inertial')
        ET.SubElement(inr, 'origin', {'xyz': vec(d['com_rel']), 'rpy': '0 0 0'})
        ET.SubElement(inr, 'mass', {'value': fmt(d['mass'])})
        ET.SubElement(inr, 'inertia', {k: fmt(d[k]) for k in ('ixx', 'ixy', 'ixz', 'iyy', 'iyz', 'izz')})
        for tag in ('visual', 'collision'):
            el = ET.SubElement(link, tag)
            ET.SubElement(el, 'origin', {'xyz': vec(d['vis']), 'rpy': vec(d['vis_rpy'])})
            g = ET.SubElement(el, 'geometry')
            s = fmt(d['scale'])
            ET.SubElement(g, 'mesh', {'filename': prefix + name + '.stl', 'scale': '%s %s %s' % (s, s, s)})
            if tag == 'visual':
                ET.SubElement(el, 'material', {'name': 'default_gray'})
    for parent, child, i in model['tree']:
        ji = jinfo[child]
        e = ji['e']
        jt = e['jt']
        if jt == 0:
            typ = 'fixed'
        elif jt == 1:
            typ = 'revolute' if ji['lo'] is not None else 'continuous'
        else:
            typ = 'prismatic'
        je = ET.SubElement(robot, 'joint', {'name': sanitize(e['name']), 'type': typ})
        ET.SubElement(je, 'parent', {'link': parent})
        ET.SubElement(je, 'child', {'link': child})
        rel = [link_world[child][k] - link_world[parent][k] for k in range(3)]
        ET.SubElement(je, 'origin', {'xyz': vec(rel), 'rpy': '0 0 0'})
        if ji['axis']:
            ET.SubElement(je, 'axis', {'xyz': vec(ji['axis'])})
        if typ in ('revolute', 'prismatic'):
            ET.SubElement(je, 'limit', {'lower': fmt(ji['lo']), 'upper': fmt(ji['hi']),
                                        'effort': fmt(DEFAULT_EFFORT), 'velocity': fmt(DEFAULT_VELOCITY)})
        elif typ == 'continuous':
            ET.SubElement(je, 'limit', {'effort': fmt(DEFAULT_EFFORT), 'velocity': fmt(DEFAULT_VELOCITY)})
    return minidom.parseString(ET.tostring(robot, 'utf-8')).toprettyxml(indent='  ')


def write_package_files(pdir, pkg):
    with open(os.path.join(pdir, 'package.xml'), 'w', encoding='utf-8') as f:
        f.write('<?xml version="1.0"?>\n<package format="3">\n  <name>%s</name>\n'
                '  <version>0.0.1</version>\n  <description>URDF exported from Fusion</description>\n'
                '  <maintainer email="user@example.com">user</maintainer>\n  <license>MIT</license>\n'
                '  <buildtool_depend>ament_cmake</buildtool_depend>\n'
                '  <export><build_type>ament_cmake</build_type></export>\n</package>\n' % pkg)
    with open(os.path.join(pdir, 'CMakeLists.txt'), 'w', encoding='utf-8') as f:
        f.write('cmake_minimum_required(VERSION 3.8)\nproject(%s)\nfind_package(ament_cmake REQUIRED)\n'
                'install(DIRECTORY urdf meshes DESTINATION share/${PROJECT_NAME})\nament_package()\n' % pkg)


# ----------------------------- STAGE 3: VERIFY -----------------------
def verify(urdf_path, pdir, rep, st):
    st['stage'], st['item'] = 'verifying the written URDF', urdf_path
    tree = ET.parse(urdf_path).getroot()
    links = {l.get('name') for l in tree.findall('link')}
    children, parents = set(), set()
    for j in tree.findall('joint'):
        p, c = j.find('parent').get('link'), j.find('child').get('link')
        if p not in links or c not in links:
            rep.error('Written joint "%s" refers to a missing link.' % j.get('name'),
                      'Report this; it is a bug in the exporter.')
        if c in children:
            rep.error('Link "%s" has two parents in the written URDF.' % c, 'Check its joints in Fusion.')
        children.add(c)
        parents.add(p)
    roots = links - children
    if len(roots) != 1:
        rep.error('The written URDF has %d root links (%s), expected 1.' % (len(roots), sorted(roots)),
                  'Check that every part is connected to base_link.')
    for l in tree.findall('link'):
        mesh = l.find('visual/geometry/mesh')
        fn = mesh.get('filename').split('/')[-1]
        p = os.path.join(pdir, 'meshes', fn)
        if not os.path.exists(p) or os.path.getsize(p) < 100:
            rep.error('Mesh file for "%s" is missing or empty.' % l.get('name'),
                      'Check that the component has visible solid bodies, then export again.')


# ----------------------------- error explanation ---------------------
def explain(e, st):
    stage, item = st.get('stage', ''), st.get('item', '')
    msg = '%s: %s' % (type(e).__name__, e)
    where = 'It failed while %s%s.' % (stage, ' (on "%s")' % item if item else '')
    fixes = {
        'reading joints': 'Delete that joint in Fusion and recreate it (Assemble > Joint), picking '
                          'geometry on two top-level components.',
        'reading joint positions and axes':
            'Edit/recreate that joint. Make sure it is a normal Joint (not As-built) and its '
            'components are not suppressed or hidden by a failed timeline item.',
        'calculating mass and inertia':
            'Check "%s" contains solid bodies (not surfaces), is not suppressed, and the timeline has '
            'no red/yellow errors.' % item,
        'exporting the STL mesh':
            'Check "%s" has visible solid bodies. Try right-click the component > Save as Mesh. If '
            'that fails too, the body needs repair (Modify > Repair Body) or simplifying.' % item,
        'creating folders':
            'Pick a different output folder you can write to (for example Documents), not a '
            'protected or synced folder.',
        'writing the URDF files':
            'Pick a different output folder, and close any program that has the old URDF open.',
        'verifying the written URDF': 'Open the written URDF in a text editor and check it is complete.',
    }
    fix = fixes.get(stage, 'Run the check again after undoing recent changes, and send me the technical detail.')
    return '%s %s (%s)' % (where, 'Fusion reported an error.', msg), fix


# ----------------------------- output --------------------------------
def finish(ui, rep, doc, headline, ok):
    text = rep.text(doc, headline)
    try:
        pal = ui.palettes.itemById('TextCommands')
        pal.isVisible = True
        pal.writeText(text)
    except:
        pass
    saved = False
    try:
        with open(REPORT_PATH, 'w', encoding='utf-8') as f:
            f.write(text)
        saved = True
    except:
        pass
    if saved and (rep.errors or rep.warnings):
        try:
            if hasattr(os, 'startfile'):
                os.startfile(REPORT_PATH)
            else:
                import subprocess
                subprocess.call(['open', REPORT_PATH])
        except:
            pass
    short = headline
    if rep.errors:
        short += '\n\nFound %d error(s) and %d warning(s).' % (len(rep.errors), len(rep.warnings))
        short += '\nEVERY error, with how to fix it, is in this file:\n' + REPORT_PATH
        short += '\n(it was opened for you; fix them all, then run again)'
    elif rep.warnings:
        short += '\n\n%d warning(s) are listed in:\n%s' % (len(rep.warnings), REPORT_PATH)
    ui.messageBox(short, 'Fusion to URDF')


def run(context):
    app = adsk.core.Application.get()
    ui = app.userInterface
    rep = Report()
    st = {'stage': 'starting', 'item': ''}
    doc = '?'
    try:
        design = adsk.fusion.Design.cast(app.activeProduct)
        if not design:
            ui.messageBox('Open a Fusion design first.', 'Fusion to URDF')
            return
        doc = design.rootComponent.name

        # 1. CHECK
        model = analyze(design, rep, st)
        if rep.errors or model is None:
            n = len(rep.errors)
            finish(ui, rep, doc, 'NOT EXPORTED: %d error(s) must be fixed first.' % n, False)
            return

        # choose folder
        dlg = ui.createFolderDialog()
        dlg.title = 'Choose where to save the URDF package'
        if dlg.showDialog() != adsk.core.DialogResults.DialogOK:
            rep.info('Export cancelled. The design passed the check.')
            finish(ui, rep, doc, 'CANCELLED: the design is ready to export.', True)
            return

        # 2. EXPORT   3. VERIFY
        res = export(design, model, dlg.folder, rep, st)
        if res is None:
            finish(ui, rep, doc, 'EXPORT FAILED on %d part(s): see the file for all of them.'
                   % len(rep.errors), False)
            return
        pdir, urdf_path, data, lw = res
        verify(urdf_path, pdir, rep, st)
        if rep.errors:
            finish(ui, rep, doc, 'EXPORTED WITH PROBLEMS: see the errors above.', False)
            return
        rep.info('Links: %d, joints: %d.' % (len(model['order']), len(model['tree'])))
        rep.info('Files: %s' % pdir)
        finish(ui, rep, doc, 'SUCCESS: URDF saved to\n%s' % urdf_path, True)

    except Exception as e:
        rep.tb = traceback.format_exc()
        p, f = explain(e, st)
        rep.error(p, f)
        finish(ui, rep, doc, 'FAILED: see the error above.', False)