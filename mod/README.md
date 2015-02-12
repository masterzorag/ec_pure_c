ec_pure_c
=========

Elliptic Curve pure C implementation

modified version to work over an OpenCL port, targeting:
- minor data requirement
- data reuse
- use of _constant, _local address spaces over _private
- precompute u8 U[20] in main and pass it as _constant