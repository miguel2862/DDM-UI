/** Fail closed on foreign model files; never reinterpret them as DDM networks. */
export function assertDDMPayload(payload: unknown): void {
  if (!payload || typeof payload !== 'object' || Array.isArray(payload)) throw new Error('Expected a DDM experiment object.');
  const data = payload as Record<string, unknown>;
  for (const key of ['model', 'modelKind']) {
    if (Array.isArray(data[key]) && data[key].length !== 1) throw new Error('Expected a single DDM model identifier.');
    const raw = Array.isArray(data[key]) ? data[key][0] : data[key];
    if (raw != null && !['dtd', 'ddm'].includes(String(raw).toLowerCase())) {
      throw new Error('This file belongs to a different model. Only DDM experiments are supported.');
    }
  }
  const layers = new Set(['US', 'PrimarySensory', 'AssociativeSensory', 'Hippocampal', 'AssociativeMotor', 'PrimaryMotor', 'Dopaminergic']);
  const checkLayer = (value: unknown) => {
    if (Array.isArray(value) && value.length !== 1) throw new Error('Expected one layer per DDM unit.');
    if (!layers.has(String(Array.isArray(value) ? value[0] : value))) throw new Error('This network contains a layer outside the DDM model.');
  };
  if (Array.isArray(data.npes)) {
    for (const value of data.npes) {
      if (!value || typeof value !== 'object') throw new Error('Expected a DDM unit object.');
      const row = value as Record<string, unknown>;
      checkLayer(row.layer ?? row.Layer);
    }
  } else if (data.npes && typeof data.npes === 'object') {
    const columns = data.npes as Record<string, unknown>;
    const column = columns.layer ?? columns.Layer;
    for (const layer of Array.isArray(column) ? column : [column]) checkLayer(layer);
  }
}
