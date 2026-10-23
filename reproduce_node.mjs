/** Reproduce the independent BigInt checks in an isolated directory. */
import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {createHash} from 'node:crypto';
import {spawnSync} from 'node:child_process';
const base=path.dirname(fileURLToPath(import.meta.url));
const source=path.join(base,'independent_node');
const args=process.argv.slice(2);
if(args.includes('--list')){console.log(['ordinary_pieri_actual_ideals.cjs','lagrangian_coefficient_windows.mjs','odd_excess_five_independent.mjs'].join('\n'));process.exit(0);}
if(args.length && (args.length!==2||args[0]!=='--output'))throw new Error('Usage: node reproduce_node.mjs [--output NEW_DIRECTORY]');
const out=path.resolve(args.length?args[1]:path.join(base,'run_node_'+new Date().toISOString().replace(/[:.]/g,'_')));
const rel=path.relative(source,out);
if(rel===''||(!rel.startsWith('..'+path.sep)&&!path.isAbsolute(rel)))throw new Error('Output must be outside the source directory');
fs.mkdirSync(out);const work=path.join(out,'work');fs.mkdirSync(work);
const names=['ordinary_pieri_actual_ideals.cjs','lagrangian_coefficient_windows.mjs','odd_excess_five_independent.mjs'];
const hash=file=>createHash('sha256').update(fs.readFileSync(file)).digest('hex');
const manifest={started_at:new Date().toISOString(),status:'RUNNING',node:process.version,scope:'Finite exact checks and a degree-bounded interpolation certificate; not formal proof-assistant verification.',source_sha256:Object.fromEntries(names.map(name=>[name,hash(path.join(source,name))])),tasks:[]};
const save=()=>fs.writeFileSync(path.join(out,'manifest.json'),JSON.stringify(manifest,null,2)+'\n');save();
try{for(const name of names)fs.copyFileSync(path.join(source,name),path.join(work,name));for(const name of names){const start=Date.now();console.log('Running '+name);const child=spawnSync(process.execPath,[path.join(work,name)],{cwd:work,windowsHide:true,encoding:'utf8',maxBuffer:8*1024*1024});fs.writeFileSync(path.join(out,name+'.log'),(child.stdout||'')+(child.stderr||''));const record={program:name,started_at:new Date(start).toISOString(),ended_at:new Date().toISOString(),exit_code:child.status,seconds:(Date.now()-start)/1000,status:child.status===0?'PASS':'FAIL'};manifest.tasks.push(record);save();if(child.error||child.status!==0)throw child.error||new Error(name+' failed');}for(const name of names)if(hash(path.join(work,name))!==manifest.source_sha256[name])throw new Error('Source changed during execution');manifest.status='PASS';manifest.ended_at=new Date().toISOString();save();console.log('All '+names.length+' Node checks passed.');}catch(error){manifest.status='FAIL';manifest.error=String(error);manifest.ended_at=new Date().toISOString();save();throw error;}
