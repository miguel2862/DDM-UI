// Test the exact TypeScript import guard without browser state or a simulator run.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const ts = require('../frontend/node_modules/typescript');
const source = fs.readFileSync(path.join(__dirname, '../frontend/src/utils/ddmCompatibility.ts'), 'utf8');
const js = ts.transpileModule(source, { compilerOptions: { module: ts.ModuleKind.CommonJS, target: ts.ScriptTarget.ES2022 } }).outputText;
const context = { exports: {} };
vm.runInNewContext(js, context);
const guard = context.exports.assertDDMPayload;
const valid = [
  {}, { model: 'dtd' }, { model: 'DDM', modelKind: 'dtd' }, { model: ['dtd'] },
  { npes: [{ layer: 'PrimaryMotor' }, { Layer: ['US'] }] },
  { npes: { Layer: ['US', 'Dopaminergic'] } },
];
const invalid = [
  null, [], { model: 'foreign' }, { modelKind: 'foreign' },
  { model: ['dtd', 'foreign'] }, { model: [] },
  { npes: [{ layer: 'Outcome' }] }, { npes: { Layer: ['US', 'Outcome'] } },
  { npes: [null] }, { npes: [{ layer: ['US', 'Outcome'] }] },
];
for (const value of valid) assert.doesNotThrow(() => guard(value));
for (const value of invalid) assert.throws(() => guard(value));
console.log('DDM import guard: ' + (valid.length + invalid.length) + ' cases passed.');
