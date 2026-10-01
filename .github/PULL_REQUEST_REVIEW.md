@copilot please review this pull request

## Context
This PR adds GuermondPopovFlux implementations for STATE and ENERGY gradient variables to extend the existing shock-capturing capability.

## Areas to focus on:
1. Mathematical correctness of the flux implementations against the Guermond-Popov (2014) formulation
2. Gradient variable transformations (especially the quotient rule for velocities in STATE case)
3. Energy correction terms consistency
4. Code style and Fortran conventions
5. Potential issues with the dependency on PR #107

Thank you for reviewing!
