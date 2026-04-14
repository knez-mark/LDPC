# LDPC
Efficient implementation of a 5G-NR LDPC encoder/decoder in C

In order to run the LDPC decoder, run the following commands:
```bash
gcc main2.c frame.c simulation.c matrix.c length_field.c LDPC_decode.c LDPC.c LDPC_encode.c LDPC_memory.c -o main
```
```bash
./main
```
