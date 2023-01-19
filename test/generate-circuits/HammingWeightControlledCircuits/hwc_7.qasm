// This circuit is a randomly generated Hamming Weight-controlled circuit. The control and target registers have 7 qubits each

OPENQASM 2.0;
include "qelib1.inc";
qreg q[23];

// Step 1: put the first register in uniform superposition, by applying hadamard to all qubits
h q[22];
h q[21];
h q[20];
h q[19];
h q[18];
h q[17];
h q[16];
// Stage 1 of Hamming control gates
// Stage 1: compute part  (7 qubits)
x q[22];
x q[21];
ccx q[22], q[21], q[15];
x q[22];
x q[21];
x q[22];
ccx q[22], q[21], q[14];
x q[22];
x q[21];
ccx q[22], q[21], q[14];
x q[21];
ccx q[22], q[21], q[13];
// Stage 2 of Hamming control gates
// Stage 2: compute part  (7 qubits)
x q[20];
ccx q[20], q[15], q[7];
x q[20];
ccx q[20], q[15], q[6];
x q[20];
ccx q[20], q[14], q[6];
x q[20];
ccx q[20], q[14], q[5];
x q[20];
ccx q[20], q[13], q[5];
x q[20];
ccx q[20], q[13], q[4];
// Stage 2: uncompute part  (n=7 qubits)
x q[22];
x q[21];
ccx q[22], q[21], q[15];
x q[22];
x q[21];
x q[22];
ccx q[22], q[21], q[14];
x q[22];
x q[21];
ccx q[22], q[21], q[14];
x q[21];
ccx q[22], q[21], q[13];
// Stage 2: swap part  (7 qubits)
//swap q[15],q[7];
cx q[15],q[7];
cx q[7],q[15];
cx q[15],q[7];
//swap q[14],q[6];
cx q[14],q[6];
cx q[6],q[14];
cx q[14],q[6];
//swap q[13],q[5];
cx q[13],q[5];
cx q[5],q[13];
cx q[13],q[5];
//swap q[12],q[4];
cx q[12],q[4];
cx q[4],q[12];
cx q[12],q[4];
// Stage 3 of Hamming control gates
// Stage 3: compute part  (7 qubits)
x q[19];
ccx q[19], q[15], q[7];
x q[19];
ccx q[19], q[15], q[6];
x q[19];
ccx q[19], q[14], q[6];
x q[19];
ccx q[19], q[14], q[5];
x q[19];
ccx q[19], q[13], q[5];
x q[19];
ccx q[19], q[13], q[4];
x q[19];
ccx q[19], q[12], q[4];
x q[19];
ccx q[19], q[12], q[3];
// Stage 3: uncompute part  (7 qubits)
Cx q[7], q[15];
x q[19];
ccx q[19], q[6], q[14];
x q[19];
ccx q[19], q[6], q[13];
x q[19];
ccx q[19], q[5], q[13];
x q[19];
ccx q[19], q[5], q[12];
x q[19];
ccx q[19], q[4], q[12];
x q[19];
ccx q[19], q[4], q[11];
ccx q[19], q[3], q[12];
// Stage 3: swap part  (7 qubits)
//swap q[15],q[7];
cx q[15],q[7];
cx q[7],q[15];
cx q[15],q[7];
//swap q[14],q[6];
cx q[14],q[6];
cx q[6],q[14];
cx q[14],q[6];
//swap q[13],q[5];
cx q[13],q[5];
cx q[5],q[13];
cx q[13],q[5];
//swap q[12],q[4];
cx q[12],q[4];
cx q[4],q[12];
cx q[12],q[4];
//swap q[11],q[3];
cx q[11],q[3];
cx q[3],q[11];
cx q[11],q[3];
// Stage 4 of Hamming control gates
// Stage 4: compute part  (7 qubits)
x q[18];
ccx q[18], q[15], q[7];
x q[18];
ccx q[18], q[15], q[6];
x q[18];
ccx q[18], q[14], q[6];
x q[18];
ccx q[18], q[14], q[5];
x q[18];
ccx q[18], q[13], q[5];
x q[18];
ccx q[18], q[13], q[4];
x q[18];
ccx q[18], q[12], q[4];
x q[18];
ccx q[18], q[12], q[3];
x q[18];
ccx q[18], q[11], q[3];
x q[18];
ccx q[18], q[11], q[2];
// Stage 4: uncompute part  (7 qubits)
Cx q[7], q[15];
x q[18];
ccx q[18], q[6], q[14];
x q[18];
ccx q[18], q[6], q[13];
x q[18];
ccx q[18], q[5], q[13];
x q[18];
ccx q[18], q[5], q[12];
x q[18];
ccx q[18], q[4], q[12];
x q[18];
ccx q[18], q[4], q[11];
x q[18];
ccx q[18], q[3], q[11];
x q[18];
ccx q[18], q[3], q[10];
ccx q[18], q[2], q[11];
// Stage 4: swap part  (7 qubits)
//swap q[15],q[7];
cx q[15],q[7];
cx q[7],q[15];
cx q[15],q[7];
//swap q[14],q[6];
cx q[14],q[6];
cx q[6],q[14];
cx q[14],q[6];
//swap q[13],q[5];
cx q[13],q[5];
cx q[5],q[13];
cx q[13],q[5];
//swap q[12],q[4];
cx q[12],q[4];
cx q[4],q[12];
cx q[12],q[4];
//swap q[11],q[3];
cx q[11],q[3];
cx q[3],q[11];
cx q[11],q[3];
//swap q[10],q[2];
cx q[10],q[2];
cx q[2],q[10];
cx q[10],q[2];
// Stage 5 of Hamming control gates
// Stage 5: compute part  (7 qubits)
x q[17];
ccx q[17], q[15], q[7];
x q[17];
ccx q[17], q[15], q[6];
x q[17];
ccx q[17], q[14], q[6];
x q[17];
ccx q[17], q[14], q[5];
x q[17];
ccx q[17], q[13], q[5];
x q[17];
ccx q[17], q[13], q[4];
x q[17];
ccx q[17], q[12], q[4];
x q[17];
ccx q[17], q[12], q[3];
x q[17];
ccx q[17], q[11], q[3];
x q[17];
ccx q[17], q[11], q[2];
x q[17];
ccx q[17], q[10], q[2];
x q[17];
ccx q[17], q[10], q[1];
// Stage 5: uncompute part  (7 qubits)
Cx q[7], q[15];
x q[17];
ccx q[17], q[6], q[14];
x q[17];
ccx q[17], q[6], q[13];
x q[17];
ccx q[17], q[5], q[13];
x q[17];
ccx q[17], q[5], q[12];
x q[17];
ccx q[17], q[4], q[12];
x q[17];
ccx q[17], q[4], q[11];
x q[17];
ccx q[17], q[3], q[11];
x q[17];
ccx q[17], q[3], q[10];
x q[17];
ccx q[17], q[2], q[10];
x q[17];
ccx q[17], q[2], q[9];
ccx q[17], q[1], q[10];
// Stage 5: swap part  (7 qubits)
//swap q[15],q[7];
cx q[15],q[7];
cx q[7],q[15];
cx q[15],q[7];
//swap q[14],q[6];
cx q[14],q[6];
cx q[6],q[14];
cx q[14],q[6];
//swap q[13],q[5];
cx q[13],q[5];
cx q[5],q[13];
cx q[13],q[5];
//swap q[12],q[4];
cx q[12],q[4];
cx q[4],q[12];
cx q[12],q[4];
//swap q[11],q[3];
cx q[11],q[3];
cx q[3],q[11];
cx q[11],q[3];
//swap q[10],q[2];
cx q[10],q[2];
cx q[2],q[10];
cx q[10],q[2];
//swap q[9],q[1];
cx q[9],q[1];
cx q[1],q[9];
cx q[9],q[1];
// Stage 6 of Hamming control gates
// Stage 6: compute part  (7 qubits)
x q[16];
ccx q[16], q[15], q[7];
x q[16];
ccx q[16], q[15], q[6];
x q[16];
ccx q[16], q[14], q[6];
x q[16];
ccx q[16], q[14], q[5];
x q[16];
ccx q[16], q[13], q[5];
x q[16];
ccx q[16], q[13], q[4];
x q[16];
ccx q[16], q[12], q[4];
x q[16];
ccx q[16], q[12], q[3];
x q[16];
ccx q[16], q[11], q[3];
x q[16];
ccx q[16], q[11], q[2];
x q[16];
ccx q[16], q[10], q[2];
x q[16];
ccx q[16], q[10], q[1];
x q[16];
ccx q[16], q[9], q[1];
x q[16];
ccx q[16], q[9], q[0];
// Stage 6: uncompute part  (7 qubits)
Cx q[7], q[15];
x q[16];
ccx q[16], q[6], q[14];
x q[16];
ccx q[16], q[6], q[13];
x q[16];
ccx q[16], q[5], q[13];
x q[16];
ccx q[16], q[5], q[12];
x q[16];
ccx q[16], q[4], q[12];
x q[16];
ccx q[16], q[4], q[11];
x q[16];
ccx q[16], q[3], q[11];
x q[16];
ccx q[16], q[3], q[10];
x q[16];
ccx q[16], q[2], q[10];
x q[16];
ccx q[16], q[2], q[9];
x q[16];
ccx q[16], q[1], q[9];
x q[16];
ccx q[16], q[1], q[8];
ccx q[16], q[0], q[9];
// Stage 6: swap part  (7 qubits)
//swap q[15],q[7];
cx q[15],q[7];
cx q[7],q[15];
cx q[15],q[7];
//swap q[14],q[6];
cx q[14],q[6];
cx q[6],q[14];
cx q[14],q[6];
//swap q[13],q[5];
cx q[13],q[5];
cx q[5],q[13];
cx q[13],q[5];
//swap q[12],q[4];
cx q[12],q[4];
cx q[4],q[12];
cx q[12],q[4];
//swap q[11],q[3];
cx q[11],q[3];
cx q[3],q[11];
cx q[11],q[3];
//swap q[10],q[2];
cx q[10],q[2];
cx q[2],q[10];
cx q[10],q[2];
//swap q[9],q[1];
cx q[9],q[1];
cx q[1],q[9];
cx q[9],q[1];
//swap q[8],q[0];
cx q[8],q[0];
cx q[0],q[8];
cx q[8],q[0];
// The Hamming-weight register has now been set up.
// Next, we add 7 Hadamard gates and 189 random controlled Clifford gates
h q[7];
h q[6];
h q[5];
h q[4];
h q[3];
h q[2];
h q[1];
// Next we add the 189 random controlled Clifford gates
cy q[13],q[6];
ccz q[10],q[5],q[2];
ccx q[12],q[2],q[6];
ccy q[10],q[3],q[7];
cs q[8],q[2];
ccz q[12],q[6],q[7];
ccz q[8],q[1],q[6];
cz q[9],q[4];
ccx q[10],q[2],q[6];
ccz q[10],q[7],q[1];
cs q[11],q[5];
ch q[13],q[6];
ccz q[14],q[2],q[5];
cx q[8],q[1];
ccy q[11],q[3],q[2];
cy q[10],q[5];
ccx q[10],q[1],q[4];
cs q[12],q[4];
cs q[13],q[1];
ccy q[8],q[3],q[7];
ccz q[15],q[1],q[3];
cy q[9],q[5];
cz q[10],q[4];
cz q[14],q[2];
ccz q[12],q[1],q[2];
cx q[11],q[4];
ccz q[9],q[2],q[6];
ch q[12],q[3];
ccy q[11],q[2],q[3];
ccx q[11],q[4],q[1];
cx q[9],q[6];
ch q[14],q[6];
cs q[14],q[6];
cs q[11],q[7];
ccy q[8],q[5],q[4];
cy q[8],q[4];
ch q[10],q[3];
cy q[15],q[6];
ccz q[13],q[5],q[4];
cy q[15],q[4];
ccx q[10],q[5],q[7];
ch q[10],q[4];
ccz q[11],q[1],q[3];
cx q[10],q[4];
ch q[15],q[1];
ccz q[8],q[2],q[5];
ccy q[15],q[2],q[4];
ch q[14],q[5];
cs q[13],q[4];
ccz q[9],q[3],q[7];
ccx q[12],q[4],q[1];
ccy q[13],q[2],q[1];
cy q[12],q[5];
ccx q[13],q[2],q[7];
ccz q[8],q[4],q[2];
cz q[8],q[2];
cz q[10],q[5];
ccy q[15],q[1],q[4];
ccx q[11],q[4],q[2];
cz q[13],q[2];
cx q[15],q[3];
cx q[11],q[3];
cx q[11],q[1];
cs q[15],q[3];
cy q[9],q[2];
cx q[14],q[4];
ccz q[9],q[3],q[4];
cy q[11],q[7];
ch q[15],q[5];
cy q[10],q[7];
cy q[9],q[3];
cx q[15],q[6];
cx q[9],q[2];
cz q[13],q[3];
cx q[14],q[1];
ccz q[11],q[3],q[7];
ccy q[15],q[4],q[7];
ccx q[8],q[3],q[5];
ccz q[10],q[4],q[3];
cy q[15],q[5];
ccx q[9],q[6],q[1];
cs q[11],q[3];
ccy q[12],q[7],q[2];
ccz q[11],q[3],q[6];
cs q[9],q[3];
cx q[15],q[4];
ccz q[15],q[2],q[7];
cx q[15],q[6];
ccx q[9],q[1],q[3];
cz q[11],q[1];
ccy q[15],q[1],q[4];
ccy q[14],q[1],q[4];
cs q[9],q[6];
cz q[10],q[3];
cx q[11],q[3];
ccz q[10],q[6],q[7];
ccz q[9],q[1],q[3];
cs q[15],q[6];
ccz q[11],q[2],q[5];
cx q[14],q[4];
cy q[10],q[6];
ccz q[10],q[4],q[1];
ccy q[10],q[3],q[5];
ccz q[15],q[5],q[3];
ccz q[8],q[4],q[6];
ccy q[8],q[7],q[3];
ccz q[13],q[6],q[1];
cs q[14],q[2];
ccx q[9],q[4],q[5];
ccy q[14],q[2],q[5];
cs q[8],q[6];
cz q[15],q[7];
ccz q[14],q[1],q[2];
cy q[8],q[2];
ccx q[11],q[7],q[5];
cx q[10],q[5];
cs q[8],q[2];
ccy q[13],q[6],q[3];
ccz q[13],q[2],q[5];
ccz q[10],q[3],q[4];
ccy q[12],q[1],q[3];
ch q[13],q[1];
ch q[12],q[1];
cz q[9],q[1];
ch q[13],q[1];
ch q[14],q[7];
ccx q[9],q[2],q[6];
cx q[15],q[4];
cx q[14],q[1];
ccx q[11],q[2],q[3];
cx q[12],q[6];
ccy q[8],q[2],q[3];
cs q[10],q[5];
cz q[12],q[6];
cs q[8],q[1];
ccz q[11],q[6],q[2];
ccy q[10],q[4],q[3];
ccx q[14],q[5],q[4];
ccx q[12],q[6],q[4];
cx q[10],q[3];
cz q[10],q[4];
ccy q[10],q[6],q[1];
cx q[10],q[3];
ccz q[12],q[4],q[1];
cs q[11],q[6];
cy q[12],q[2];
ccy q[12],q[5],q[2];
cs q[9],q[5];
ccx q[12],q[2],q[5];
ccz q[8],q[1],q[5];
ccz q[12],q[3],q[6];
ccx q[9],q[1],q[2];
cz q[9],q[2];
ccy q[13],q[7],q[6];
ch q[10],q[1];
ccx q[15],q[7],q[3];
ccz q[8],q[7],q[2];
ccx q[13],q[7],q[2];
ch q[10],q[6];
ch q[8],q[2];
cx q[11],q[4];
ccx q[12],q[5],q[2];
cx q[13],q[6];
cs q[10],q[4];
ccz q[13],q[1],q[4];
cx q[8],q[1];
ccx q[11],q[5],q[7];
ccx q[8],q[3],q[6];
cy q[14],q[2];
ccy q[15],q[5],q[6];
ccz q[8],q[2],q[1];
cz q[14],q[1];
ccy q[11],q[2],q[1];
ccy q[10],q[3],q[2];
cy q[14],q[7];
ccy q[11],q[7],q[5];
cy q[14],q[7];
cy q[13],q[1];
cz q[10],q[1];
cx q[12],q[7];
ccy q[12],q[7],q[2];
ccy q[13],q[5],q[3];
cs q[8],q[6];
cx q[8],q[3];
cy q[14],q[2];
cz q[14],q[6];
ccx q[10],q[4],q[6];
ccy q[9],q[3],q[2];
ch q[13],q[6];
