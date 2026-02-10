# synNotch-Gillespie-Contact-Variation
Outlines and codes for Gillespie simulations of synNotch synthesis with time-varied contact

## The Gillespie Algorithm
### Third Level Heading



$$
\frac{dG_\alpha}{dt} = \frac{l_{\alpha \beta}}{P_{\alpha}} - D G_\alpha
$$


### Reaction 1: Production
$\emptyset \xrightarrow{k_{\text{prod}}} G$

$$
\emptyset \xrightarrow{k_{\text{prod}}} G
$$

with rate:
$$
k_{\text{prod}} = \frac{S L}{P}
$$

### Reaction 2: Degradation
$$
G \xrightarrow{D} \emptyset
$$

with rate:
$$
k_{\text{deg}} = D \cdot G
$$


# References
http://www.be150.caltech.edu/2019/handouts/12_stochastic_simulation_all_code.html
