const ZERO=Object.freeze({re:0,im:0});
export class ACAnalysisError extends Error{constructor(message){super(message);this.name='ACAnalysisError';}}

export function analyzeAC(circuit,options={}){
  validateCircuit(circuit);
  const frequencies=createFrequencySweep(options);
  const ground=circuit.ground??'0';
  const used=new Set(circuit.components.flatMap(component=>[component.from,component.to]));
  const nodes=[...used].filter(node=>node!==ground);
  const nodeIndex=new Map(nodes.map((node,index)=>[node,index]));
  const voltageSources=circuit.components.filter(component=>component.type==='voltage-source');
  const sourceIndex=new Map(voltageSources.map((source,index)=>[source.id,nodes.length+index]));
  const size=nodes.length+voltageSources.length;
  const result={analysis:'ac',scale:options.scale??'log',frequencies,nodeVoltages:Object.fromEntries([ground,...nodes].map(node=>[node,[]])),componentCurrents:Object.fromEntries(circuit.components.map(component=>[component.id,[]]))};
  const ni=node=>node===ground?undefined:nodeIndex.get(node);
  for(const frequency of frequencies){
    const omega=2*Math.PI*frequency;
    const matrix=Array.from({length:size},()=>Array.from({length:size},()=>complex()));
    const rhs=Array.from({length:size},()=>complex());
    const stampY=(from,to,y)=>{const i=ni(from),j=ni(to);if(i!==undefined)matrix[i][i]=add(matrix[i][i],y);if(j!==undefined)matrix[j][j]=add(matrix[j][j],y);if(i!==undefined&&j!==undefined){matrix[i][j]=sub(matrix[i][j],y);matrix[j][i]=sub(matrix[j][i],y);}};
    const stampI=(from,to,current)=>{const i=ni(from),j=ni(to);if(i!==undefined)rhs[i]=sub(rhs[i],current);if(j!==undefined)rhs[j]=add(rhs[j],current);};
    for(const component of circuit.components){
      if(component.type==='resistor')stampY(component.from,component.to,complex(1/component.resistance));
      else if(component.type==='capacitor')stampY(component.from,component.to,complex(0,omega*component.capacitance));
      else if(component.type==='inductor')stampY(component.from,component.to,complex(0,-1/(omega*component.inductance)));
      else if(component.type==='current-source')stampI(component.from,component.to,phasor(component));
      else if(component.type==='voltage-source'){
        const k=sourceIndex.get(component.id),i=ni(component.from),j=ni(component.to);
        if(i!==undefined){matrix[i][k]=add(matrix[i][k],complex(1));matrix[k][i]=add(matrix[k][i],complex(1));}
        if(j!==undefined){matrix[j][k]=sub(matrix[j][k],complex(1));matrix[k][j]=sub(matrix[k][j],complex(1));}
        rhs[k]=phasor(component);
      }
    }
    const solution=solveComplex(matrix,rhs);
    const voltage=node=>node===ground?ZERO:solution[nodeIndex.get(node)];
    for(const node of [ground,...nodes])result.nodeVoltages[node].push(clean(voltage(node)));
    for(const component of circuit.components){
      const v=sub(voltage(component.from),voltage(component.to));let current;
      if(component.type==='resistor')current=scale(v,1/component.resistance);
      else if(component.type==='capacitor')current=mul(v,complex(0,omega*component.capacitance));
      else if(component.type==='inductor')current=mul(v,complex(0,-1/(omega*component.inductance)));
      else if(component.type==='current-source')current=phasor(component);
      else current=solution[sourceIndex.get(component.id)];
      result.componentCurrents[component.id].push(clean(current));
    }
  }
  return result;
}

export function createFrequencySweep({startFrequency=10,stopFrequency=1e6,points=101,scale:kind='log'}={}){
  if(!Number.isFinite(startFrequency)||startFrequency<=0)throw new ACAnalysisError('startFrequency must be positive');
  if(!Number.isFinite(stopFrequency)||stopFrequency<startFrequency)throw new ACAnalysisError('stopFrequency must be at least startFrequency');
  if(!Number.isInteger(points)||points<1||points>100000)throw new ACAnalysisError('points must be an integer from 1 to 100000');
  if(!['log','linear'].includes(kind))throw new ACAnalysisError('scale must be "log" or "linear"');
  if(points===1)return[startFrequency];
  return Array.from({length:points},(_,index)=>kind==='linear'?startFrequency+(stopFrequency-startFrequency)*index/(points-1):startFrequency*Math.pow(stopFrequency/startFrequency,index/(points-1)));
}

export function toBode(series){return{magnitude:series.map(abs),magnitudeDb:series.map(value=>20*Math.log10(abs(value))),phaseDegrees:series.map(value=>Math.atan2(value.im,value.re)*180/Math.PI)};}

export function transferFunction(result,outputNode,inputNode){const output=result.nodeVoltages[outputNode],input=result.nodeVoltages[inputNode];if(!output||!input)throw new ACAnalysisError('outputNode and inputNode must exist in the result');return output.map((value,index)=>div(value,input[index]));}

export function findCutoffFrequency(frequencies,series,{dropDb=3,reference='maximum'}={}){
  if(frequencies.length!==series.length||frequencies.length<2)throw new ACAnalysisError('frequency and series arrays must have matching lengths of at least two');
  const db=toBode(series).magnitudeDb;const base=reference==='first'?db[0]:Math.max(...db);const target=base-dropDb;
  for(let index=1;index<db.length;index++)if(db[index]<=target&&db[index-1]>target){const ratio=(target-db[index-1])/(db[index]-db[index-1]);const logFrequency=Math.log(frequencies[index-1])+ratio*(Math.log(frequencies[index])-Math.log(frequencies[index-1]));return Math.exp(logFrequency);}
  return undefined;
}

export function solveComplex(matrix,vector,tolerance=1e-14){
  const n=vector.length;if(n===0)return[];if(matrix.length!==n||matrix.some(row=>row.length!==n))throw new ACAnalysisError('matrix dimensions must match vector');
  const augmented=matrix.map((row,index)=>[...row.map(copy),copy(vector[index])]);
  for(let column=0;column<n;column++){
    let pivot=column;for(let row=column+1;row<n;row++)if(abs(augmented[row][column])>abs(augmented[pivot][column]))pivot=row;
    if(abs(augmented[pivot][column])<=tolerance)throw new ACAnalysisError('circuit has no unique AC solution; check floating nodes and source conflicts');
    [augmented[column],augmented[pivot]]=[augmented[pivot],augmented[column]];
    const divisor=augmented[column][column];for(let index=column;index<=n;index++)augmented[column][index]=div(augmented[column][index],divisor);
    for(let row=0;row<n;row++){if(row===column)continue;const factor=augmented[row][column];for(let index=column;index<=n;index++)augmented[row][index]=sub(augmented[row][index],mul(factor,augmented[column][index]));}
  }
  return augmented.map(row=>clean(row[n]));
}

export function complex(re=0,im=0){return{re,im};}
export function add(a,b){return complex(a.re+b.re,a.im+b.im);}
export function sub(a,b){return complex(a.re-b.re,a.im-b.im);}
export function mul(a,b){return complex(a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re);}
export function div(a,b){const denominator=b.re*b.re+b.im*b.im;if(denominator===0)throw new ACAnalysisError('division by zero phasor');return complex((a.re*b.re+a.im*b.im)/denominator,(a.im*b.re-a.re*b.im)/denominator);}
export function scale(a,factor){return complex(a.re*factor,a.im*factor);}
export function abs(a){return Math.hypot(a.re,a.im);}

function phasor(component){const magnitude=component.acMagnitude??0,phase=(component.acPhase??0)*Math.PI/180;return complex(magnitude*Math.cos(phase),magnitude*Math.sin(phase));}
function copy(value){return complex(value.re,value.im);}
function clean(value){return complex(Math.abs(value.re)<1e-15?0:value.re,Math.abs(value.im)<1e-15?0:value.im);}
function validateCircuit(circuit){
  if(!circuit||!Array.isArray(circuit.nodes)||!Array.isArray(circuit.components))throw new ACAnalysisError('circuit requires nodes and components arrays');
  const ids=new Set(circuit.nodes.map(node=>node.id)),ground=circuit.ground??'0';if(!ids.has(ground))throw new ACAnalysisError('ground node is missing');
  const supported=new Set(['resistor','capacitor','inductor','voltage-source','current-source']);const componentIds=new Set();
  for(const component of circuit.components){if(!component.id||componentIds.has(component.id))throw new ACAnalysisError('component ids must be unique');componentIds.add(component.id);if(!supported.has(component.type))throw new ACAnalysisError(`unsupported component type "${component.type}"`);if(!ids.has(component.from)||!ids.has(component.to))throw new ACAnalysisError(`${component.id} references an unknown node`);if(component.from===component.to)throw new ACAnalysisError(`${component.id} connects a node to itself`);for(const field of ['resistance','capacitance','inductance'])if(component[field]!==undefined&&(!Number.isFinite(component[field])||component[field]<=0))throw new ACAnalysisError(`${component.id}.${field} must be positive`);if(component.acMagnitude!==undefined&&(!Number.isFinite(component.acMagnitude)||component.acMagnitude<0))throw new ACAnalysisError(`${component.id}.acMagnitude must be non-negative`);}
}
