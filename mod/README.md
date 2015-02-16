ec_pure_c
=========

Elliptic Curve pure C implementation

modified version to work over an OpenCL port, targeting:
- minor data requirement
- data reuse
- use of _constant, _local address spaces over _private
- precompute U, V needed in bn_mon_inv, move to _constant
