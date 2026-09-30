import {abs,analyzeAC,findCutoffFrequency,transferFunction} from '../src/index.js';
const circuit={ground:'0',nodes:[{id:'0'},{id:'in'},{id:'out'}],components:[{id:'V1',type:'voltage-source',from:'in',to:'0',acMagnitude:1},{id:'R1',type:'resistor',from:'in',to:'out',resistance:1000},{id:'C1',type:'capacitor',from:'out',to:'0',capacitance:159.154943e-9}]};
const result=analyzeAC(circuit,{startFrequency:10,stopFrequency:100000,points:401});
const transfer=transferFunction(result,'out','in');
console.log(JSON.stringify({cutoffFrequency:findCutoffFrequency(result.frequencies,transfer),gainAt10Hz:abs(transfer[0]),gainAt100kHz:abs(transfer.at(-1))},null,2));
