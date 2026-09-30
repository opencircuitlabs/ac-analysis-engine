# OpenCircuit AC Analysis Engine

A dependency-free small-signal AC circuit simulator using complex modified nodal analysis. It produces frequency-domain node voltages and component currents suitable for Bode plots, impedance analysis, and filter verification.

## Capabilities

- Resistors, capacitors, inductors, independent voltage sources, and independent current sources.
- Linear or logarithmic sweeps of up to 100,000 points.
- Complex Gaussian elimination with partial pivoting.
- Source magnitude and phase excitation.
- Node-voltage and branch-current phasors at every frequency.
- Bode magnitude, decibel, and phase conversion.
- Transfer-function and interpolated cutoff-frequency helpers.
- Explicit validation and singular-circuit diagnostics.
- No native modules or runtime dependencies.

## Installation

```bash
npm install github:opencircuitlabs/ac-analysis-engine
```

## Example

```js
import { analyzeAC, toBode, transferFunction } from '@opencircuitlabs/ac-analysis-engine';

const circuit = {
  ground: '0',
  nodes: [{id:'0'}, {id:'in'}, {id:'out'}],
  components: [
    {id:'V1', type:'voltage-source', from:'in', to:'0', acMagnitude:1},
    {id:'R1', type:'resistor', from:'in', to:'out', resistance:1000},
    {id:'C1', type:'capacitor', from:'out', to:'0', capacitance:159.154943e-9}
  ]
};

const result = analyzeAC(circuit, {
  startFrequency: 10,
  stopFrequency: 100000,
  points: 401,
  scale: 'log'
});

const bode = toBode(transferFunction(result, 'out', 'in'));
```

Phasors use portable objects in the form `{ re, im }`. Component values use SI units and source phase is expressed in degrees.

## Numerical method

At each frequency, the engine stamps complex admittances into a modified nodal analysis matrix. Voltage sources introduce branch-current unknowns. The matrix is solved with complex Gaussian elimination and magnitude-based partial pivoting.

The engine performs linear small-signal analysis. Nonlinear semiconductor devices must be linearized around a DC operating point by a future operating-point engine before being passed here.

## Development

```bash
npm test
npm run example
```

Tests cover analytic divider, RC cutoff, phase, transfer, sweep, pivoting, singularity, and validation behavior.

## License

MIT © 2026 OpenCircuitLabs
