import { runAudit } from './two_row_symbolic_referee.mjs';
import { readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
const source = process.argv[2];
const result = runAudit();
if (!result.pass) throw new Error('Exact symbolic identity failed');
const hash = async (url) => createHash('sha256').update(await readFile(url)).digest('hex');
const record = {
  sourceSha256: source ? await hash(source) : null,
  checkerSha256: await hash(new URL('./two_row_symbolic_referee.mjs', import.meta.url)),
  ...result
};
await writeFile(new URL('./two_row_symbolic_referee.json', import.meta.url), JSON.stringify(record, null, 2)+'\n');
console.log(JSON.stringify(record, null, 2));
