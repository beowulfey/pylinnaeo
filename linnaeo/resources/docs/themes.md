## THEMES

Linnaeo provides a variety of different themes for displaying amino acid properties in the sequence viewer. 
This document is designed to provide a deeper explanation for the meaning behind what these color schemes represent.

---
### Linnaeo default theme

The default theme was partially inspired by Clustal X categories but with colors that were intended to be as easily readable as possible. Although not necessarily the best for quickly seeing unusual residues or differences, it may perhaps be the easiest on the eyes. 

Colors used for the default theme: 

<span style="background-color: #97a4e8; color: #000000">Hydrophobic</span>, 
<span style="background-color: #a0edd8; color: #000000">Aromatic</span>, 
<span style="background-color: #bef1ac; color: #000000">Polar</span>, 
<span style="background-color: #f4f2ba; color: #000000">Cysteine</span>, 
<span style="background-color: #f6decc; color: #000000">Proline</span>, 
<span style="background-color: #db8a8b; color: #000000">Positive</span>, 
<span style="background-color: #e190e2; color: #000000">Negative</span>, 
<span style="background-color: #f7edec; color: #000000">Alanine</span>, 
<span style="background-color: #efefef; color: #000000">Glycine</span>

Residues are organized into the following categories, based on property: 

|Hydrophobic|Aromatic|Polar|Positive|Negative|Misc|
|---|---|---|---|---|---|
|Leucine (L)|Tryptophan (W)|Serine (S)|Arginine (R)|Aspartate (D)|Cysteine (C)|
|Isoleucine (I)|Tyrosine (Y)|Threonine (T)|Lysine (K)|Glutamate (E)|Proline (P)
|Valine (V)|Phenylalanine (F)| Aspargine (N)|Histidine (H)||Alanine (A)|
|Methionine (M)||Glutamine (Q)|||Glycine (G)|

This same organization is used by the following color schemes: 
* Bold theme (<span style=" background-color: #6c7edf; color: #000000">Hydrophobic</span>, 
<span style=" background-color: #79eaca; color: #000000">Aromatic</span>, 
<span style=" background-color: #a5f28b; color: #000000">Polar</span>, 
<span style=" background-color: #f7f5a1; color: #000000">Cysteine</span>, 
<span style=" background-color: #fbd7bc; color: #000000">Proline</span>, 
<span style=" background-color: #c4585a; color: #000000">Positive</span>, 
<span style=" background-color: #d161d2; color: #000000">Negative</span>, 
<span style=" background-color: #feefee; color: #000000">Alanine</span>, 
<span style=" background-color: #efefef; color: #000000">Glycine</span>)
* Grayscale (<span style=" background-color: #d6d6d6; color: #000000">Hydrophobic</span>,
<span style=" background-color: #000000; color: #FFFFFF">Aromatic</span>,
<span style=" background-color: #b3b3b3; color: #000000">Polar</span>,
<span style=" background-color: #242424; color: #FFFFFF">Cysteine</span>,
<span style=" background-color: #474747; color: #FFFFFF">Proline</span>,
<span style=" background-color: #8f8f8f; color: #000000">Positive</span>,
<span style=" background-color: #6b6b6b; color: #FFFFFF">Negative</span>,
<span style=" background-color: #e8e8e8; color: #000000">Alanine</span>,
<span style=" background-color: #fafafa; color: #000000">Glycine</span>)
* Rainbow (<span style=" background-color: #d0ff22; color: #000000">Hydrophobic</span>,
<span style=" background-color: #22ffc0; color: #000000">Aromatic</span>,
<span style=" background-color: #eae9e9; color: #000000">Alanine</span>,
<span style=" background-color: #22cfff; color: #000000">Negative</span>,
<span style=" background-color: #ff5122; color: #000000">Positive</span>,
<span style=" background-color: #ffbf22; color: #000000">Polar</span>,
<span style=" background-color: #c022ff; color: #000000">Cysteine</span>,
<span style=" background-color: #62ff22; color: #000000">Proline</span>,
<span style=" background-color: #eae9e9; color: #000000">Glycine</span>)
* Colorsafe (<span style=" background-color: #a1c8ec; color: #000000">Hydrophobic</span>,
<span style=" background-color: #5f9ed1; color: #000000">Aromatic</span>,
<span style=" background-color: #0167ab; color: #FFFFFF">Alanine</span>,
<span style=" background-color: #fbbd77; color: #000000">Negative</span>,
<span style=" background-color: #ff800f; color: #000000">Positive</span>,
<span style=" background-color: #c85200; color: #FFFFFF">Polar</span>,
<span style=" background-color: #898989; color: #000000">Cysteine</span>,
<span style=" background-color: #cfcfcf; color: #000000">Proline</span>,
<span style=" background-color: #efefef; color: #000000">Glycine</span>)

---

### Hydropathy theme

This theme colors each individual amino acid based on its calculated preference for water by the Kyle-Doolittle scale. 
More information on this is available [here](http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html#:~:text=The%20Kyte%2DDoolittle%20scale%20is,on%20the%20window%20size%20used.) 
and [here](http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html).

The colors range in a continuous scale from blue to red, via purple, where blue is the most hydrophobic and red is the most hydrophilic. 
Note that although the color scale varies, the amount each color varies by does not perfectly correlate with the changing amino acid values on the K-D scale. 
[See the tool used to generate the colors here](https://bit.ly/2V8VF11).

||||||||
|---|---|---|---|---|---|---|
|<span style=" background-color: #6FA6FF; color: #000000">&nbsp;Isoleucine&nbsp;</span>|<span style=" background-color: #6FA5FF; color: #000000">&nbsp;Valine&nbsp;</span>|<span style=" background-color: #6FA4FF; color: #000000">&nbsp;Leucine&nbsp;</span>|<span style=" background-color: #6F9EFF; color: #000000">&nbsp;Phenylalanine&nbsp;</span>|<span style=" background-color: #6F94FF; color: #000000">&nbsp;Cysteine&nbsp;</span>|<span style=" background-color: #6E86FE; color: #000000">&nbsp;Methionine&nbsp;</span>|<span style=" background-color: #6E74FE; color: #000000">&nbsp;Alanine&nbsp;</span>|
|<span style=" background-color: #7D6DFD; color: #000000">&nbsp;Tryptophan&nbsp;</span>|<span style=" background-color: #966CFC; color: #000000">&nbsp;Glycine&nbsp;</span>|<span style=" background-color: #B26BFA; color: #000000">&nbsp;Threonine&nbsp;</span>|<span style=" background-color: #CF69F6; color: #000000">&nbsp;Serine&nbsp;</span>|<span style=" background-color: #DB60E2; color: #000000">&nbsp;Tyrosine&nbsp;</span>|<span style=" background-color: #DE5EC9; color: #000000">&nbsp;Proline&nbsp;</span>|<span style=" background-color: #DC5DAF; color: #000000">&nbsp;Histidine&nbsp;</span>|||
|<span style=" background-color: #DA5C86; color: #000000">&nbsp;Asparagine&nbsp;</span>|<span style=" background-color: #DA5C78; color: #000000">&nbsp;Aspartate&nbsp;</span>|<span style=" background-color: #D95B6D; color: #000000">&nbsp;Glutamine&nbsp;</span>|<span style=" background-color: #D95B65; color: #000000">&nbsp;Glutamate&nbsp;</span>|<span style=" background-color: #D95B5F; color: #000000">&nbsp;Lysine&nbsp;</span>|<span style=" background-color: #D95A5D; color: #000000">&nbsp;Arginine&nbsp;</span>|

---

### Conservation theme

This theme calculates the relative similarity between a residue and the equivalent residue of the chosen reference sequence. The determination results in one of three colors depending on the resulting category: 
<span style=" background-color: #ffc380; color: #000000">Same property</span>,
<span style=" background-color: #97d3aa; color: #000000">Compatible</span>,
<span style=" background-color: #96c3d6; color: #000000">Dissimilar</span>, or
<span style=" color: #000000">Not conserved</span>

If a RESIDUE and the REFERENCE residue are both contained within one of the following sets, they fall into the respective category. The sets are searched in order, so if they aren't found in the first category, it will then search the second, or third, etc. 
These categories are constructed based on a number of different criteria. Some of these criteria can be [found here](http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html).


**<span style=" background-color: #ffc380; color: #000000">&nbsp;Same property:&nbsp;</span>**

These changes are very similar in property, and likely to have little effect on the protein as a whole. 
* Residues are identical
* {A, I, L, V, M}
* {R, H, K}
* {S, T}
* {D, E}
* {N, Q}

**<span style=" background-color: #97d3aa; color: #000000">&nbsp;Compatible:&nbsp;</span>**

These changes are mostly compatible, and may require compensatory changes elsewhere in the protein but should still be similar overall 
* {I, L, V, M, F}
* {F, W, Y}
* {A, S}
* {D, N}
* {E, Q}
* {A, T}
* {H, Q}

**<span style=" background-color: #96c3d6; color: #000000">&nbsp;Dissimilar:&nbsp;</span>**

These changes are suggestive of larger functional changes and less overall conservation. 
* {T, I, L, V, M, F}
* {R, E, Q}
* {K, E, Q}
* {C, S}
* {H, Y}
* {W, Q}

**<span style=" color: #000000">&nbsp;Not conserved:&nbsp;</span>**

Any changes not contained within the above sets are considered functionally unconserved and are not colored.


Note that this is providing a sense of conservation based on residue *properties*... to view conservation based on the exact sequence itself, 
the "Use reference" option can be turned on for any other theme. With this setting, a residue will remain uncolored if it differs from the 
reference amino acid at all. 
