[![DOI](https://zenodo.org/badge/492541524.svg)](https://zenodo.org/doi/10.5281/zenodo.10080462)


This repository contains the accompanying codes for paper [Optimizing implementations of linear layers using two and higher input XOR gates](https://peerj.com/articles/cs-1820/#).

**Cite this as:**
```

Kurt Pehlivanoğlu M, Demir MA. 2024. Optimizing implementations of linear layers using two and higher input XOR gates. PeerJ Computer Science 10:e1820 https://doi.org/10.7717/peerj-cs.1820

```



# SBP (Superior Boyar-Peralta) 
# dB-BDKCI (depth-bounded version of BDKCI)

SBP heuristic algorithm provides global optimization solutions, especially for extracting low-latency circuits.

The dB-BDKCI heuristic algorithm improves upon the [BDKCI heuristic algorithm](https://eprint.iacr.org/2021/1400), which was recently proposed, by integrating circuit depth awareness and constraining the depth of the circuits. Using the dB-BDKCI heuristic, we present better circuit implementations of linear layers of block ciphers than those given in the literature. Moreover, the given circuit for the AES MixColumn matrix only requires 44 XOR gates/depth 3/240.95 GE in STM 130nm library.
