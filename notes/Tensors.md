---
created: "2025-12-17T15:24:01.287251"
modified: "2026-07-07T22:58:14.345555"
title: "Tensors"
---

# Tensors

### Scalars

In physics, a _scalar_ is a quantity that can be fully described by just one measurable number, with associated physical units. Examples include temperature, mass, and volume. We represent scalars with italic-face symbols, like $T$ for temperature, $m$ for mass, and $V$ for volume.

The "one number" definition is meant to contrast scalars with the higher-order objects we'll discuss next. Scalars can vary with space (scalar _fields_) and/or time, and the measurement we ascribe to them can depend on scale. Temperature for example is a macroscopic quantity that starts to lose meaning at sufficiently microscopic scales. When we later discuss differential calculus of scalars, approximation to idealized _continuous_ scalar fields will prove convenient.

### Vectors

A _vector_ is a quantity that has direction and magnitude (_i.e._, length). We write vectors as roman-face symbols like $\mathrm{r}$ for position, $\mathrm{v}$ for velocity, and $\mathrm{E}$ for an electric field. We can visualize vectors as arrows in our three dimensions.

The magnitude of a vector does not depend on, or is _invariant_ to, the way we frame the space in which that the vector resides, _i.e._, the choice of the spatial coordinates. A vector's magnitude is a scalar and we can write it like

$$|\mathrm{v}|=v$$

which says the velocity vector $\mathrm{v}$ has a magnitude equal to a scalar $v$ known as the _speed_.

The purely directional part of a vector is the _normalized_ (length of 1) unit vector we denote with double brackets:

$$\llbracket \mathrm{v} \rrbracket = \frac{\mathrm{v}}{|\mathrm{v}|}$$

Unlike its magnitude, a vector's measured _direction_ does depend on the choice of spatial coordinates. Having already specified the vector's magnitude, a measurement of direction in our three-dimensional space can consist of up to two numbers, _e.g._, the direction of the position vector $\mathrm{r}$ can be measured as angular components $\phi, \theta$ when we use a spherical coordinate system. The actual measured values of $\phi$ and $\theta$ that we'd obtain depend on how the coordinate system is defined and oriented.

In addition to spherical coordinates, another typical coordinate system of choice is of the $x,y,z$ Cartesian variety. We can write

$$\mathrm{r}=x\mathrm{i}+y\mathrm{j}+z\mathrm{k}$$

where $x,y,z$ are scalar components and $\mathrm{i},\mathrm{j},\mathrm{k}$ are orthogonal unit vectors (_unit_ meaning magnitude of 1) comprising the coordinate basis. Instead of these particular symbols we can write

$$\mathrm{r}=\sum_{i=1}^{3}{x_i\mathrm{e}_i}$$

or in the more compact shorthand:

$$\mathrm{r}=x_i\mathrm{e}_i$$

or the even more compact shorthand:

$$\mathrm{r}=x_i$$

The single index $i$ is understood in this more compact notation to run over our three spatial dimensions and include three scalar components each multiplied by a unit basis vector. This last notational format was established by Albert Einstein in 1916 and is therefore known as "Einstein notation".

Irrespective of how we use indices in the notation &mdash; or indeed, whether we choose to use indices at all &mdash; a vector can still be thought of as having an intrinsic direction that exists independent of whatever coordinates we choose to measure it against. A primary motivation of this work is to arrive at, and appreciate, index-free coordinate-free notations. It's also undeniable that simply writing $\mathrm{r}$ is by far the _most_ compact representation out of all of the above (so compact that we really have to pay attention to the roman typeface as the sole notational indication that it's a vector)!

Addition and subtraction of vectors also has a coordinate-free meaning and can be visualized geometrically without the need for any coordinate frame.

### Dyads

Considering two vectors $\mathrm{v}$ and $\mathrm{w}$ as a pair, they define an _outer product_
$$\mathsf{D}=\mathrm{v}\mathrm{w}$$
where $\mathsf{D}$ is referred to as a _second-order tensor_, in this case specifically a _dyad_. We use san-serif typeface to represent second-order tensors in an index-free manner.

For now, we'll hold off on discussing what this dyad object actually represents and how we can or should interpret it from a physical standpoint. Instead, for the moment we'll examine the mathematical structure it has, using the same coordinate system scaffolding that applies to vectors.

###  Dyad coordinate bases

To understand how the dyad $$\mathsf{D}=\mathrm{v}\mathrm{w}$$ is structured, let's look at it in the Cartesian space:

$$\mathsf{D}=(v_x\mathrm{i}+v_y\mathrm{j}+v_z\mathrm{k})(w_x\mathrm{i}+w_y\mathrm{j}+w_z\mathrm{k})\\=v_x w_x\mathrm{ii}+v_x w_y \mathrm{ij}+v_x w_z\mathrm{ik}+v_y w_x\mathrm{ji}+v_y w_y \mathrm{jj}+...$$

There are nine total basis dyads (the outer product combinations of $\mathrm{i},\mathrm{j},\mathrm{k}$) accompanied by nine total scalar components. In index form,

$$D_{ij}=v_i w_j$$

We can also represent it as a $3\times 3$ matrix:

$$\mathsf{D}=\left( \begin{array}{ccc} v_x w_x & v_x w_y & v_x w_z \\ v_y w_x & v_y w_y & v_y w_z \\ v_z w_x & v_z w_y & v_z w_z \end{array} \right)$$

This matrix form is still interpretable as the sum of nine basis matrices, related to the above basis dyads according to:

$$
\begin{aligned}
\mathrm{ii} &=
\begin{pmatrix}
 1 & 0 & 0\\
 0 & 0 & 0\\
 0 & 0 & 0
\end{pmatrix}\quad
\mathrm{ij} =
\begin{pmatrix}
 0 & 1 & 0\\
 0 & 0 & 0\\
 0 & 0 & 0
\end{pmatrix}
\end{aligned}
$$

and so on for the others.

Given a general matrix 

$$\mathsf{A}=\begin{pmatrix}
 a & b & c\\
 d & e & f\\
 g & h & i
\end{pmatrix},$$

its _transpose_ is defined as the matrix whose elements are "reflected" about the diagonal:

$$\mathsf{A}^\mathsf{T}=\begin{pmatrix}
 a & d & g\\
 b & e & h\\
 c & f & i
\end{pmatrix}.$$

The _symmetric part_ of the matrix is:

$$\mathsf{A}_{\mathrm{symm}}=\frac{1}{2}(\mathsf{A}+\mathsf{A}^\mathsf{T})$$
$$=\begin{pmatrix}
 a & \frac{1}{2}(b+d) & \frac{1}{2}(c+g)\\
  \frac{1}{2}(d+b) & e & \frac{1}{2}(f+h)\\
 \frac{1}{2}(g+c)& \frac{1}{2}(h+f) & i
\end{pmatrix}$$

which has 6 unique components, while the _antisymmetric part_ is:

$$\mathsf{A}_{\mathrm{anti}}=\frac{1}{2}(\mathsf{A}-\mathsf{A}^\mathsf{T})$$
$$=\begin{pmatrix}
 0 & \frac{1}{2}(b-d) & \frac{1}{2}(c-g)\\
  \frac{1}{2}(d-b) & 0 & \frac{1}{2}(f-h)\\
 \frac{1}{2}(g-c)& \frac{1}{2}(h-f) & 0
\end{pmatrix}$$

which has 3 unique components, not counting the sign differences on other side of the diagonal.

Lastly, note that $\mathsf{A}$ is the sum of the symmetric and antisymmetric parts (which is why we're calling them "parts" in the first place).

For the specific case of our dyad $\mathsf{D}=\mathrm{vw}$, the transpose is such that $\mathsf{D}^\mathsf{T}=\mathrm{wv}$, and 

$$\mathsf{D}_\mathrm{symm}=\frac{1}{2}(\mathrm{vw}+\mathrm{wv})$$

$$\mathsf{D}_\mathrm{anti}=\frac{1}{2}(\mathrm{vw}-\mathrm{wv})$$

All of this guides us to a particular set of basis dyads that cleanly separate six symmetric components from three antisymmetric components. The choice of symmetric basis dyads, based on the matrix representation, is

$$\begin{array}{ccc} \mathrm{ii}, & \mathrm{jj}, & \mathrm{kk}, \\ \frac{1}{\sqrt{2}}(\mathrm{ij} +\mathrm{ji}), & \frac{1}{\sqrt{2}}(\mathrm{jk} + \mathrm{kj}), & \frac{1}{\sqrt{2}}(\mathrm{ki} +\mathrm{ik}). \end{array}$$

The remaining three basis dyads, for the antisymmetric part, are

$$\begin{array}{ccc} \frac{1}{\sqrt{2}}(\mathrm{ij} -\mathrm{ji}), & \frac{1}{\sqrt{2}}(\mathrm{jk} - \mathrm{kj}), & \frac{1}{\sqrt{2}}(\mathrm{ki} - \mathrm{ik}). \end{array}$$

The six components of $\mathsf{D}$ for the symmetric part are

$$\begin{array}{cc} v_x w_x, & \frac{1}{\sqrt{2}}(v_x w_y + v_y w_x), \\ v_y w_y, & \frac{1}{\sqrt{2}}(v_y w_z + v_z w_y), \\ v_z w_z, & \frac{1}{\sqrt{2}}(v_z w_x +v_x w_z) \end{array}$$

and the antisymmetric components are 

$$\begin{array}{c} \frac{1}{\sqrt{2}}(v_x w_y - v_y w_x), \\ \frac{1}{\sqrt{2}}(v_y w_z - v_z w_y), \\ \frac{1}{\sqrt{2}}(v_z w_x - v_x w_z). \end{array}$$

With some creative manipulation of terms, we can replace the contribution from the first three basis elements:

$$v_x w_x \mathrm{ii}+v_y w_y \mathrm{jj}+v_z w_z\mathrm{kk},$$

with the equivalent expression 

$$\frac{1}{3}(v_x w_x+v_y w_y+v_z w_z)(\mathrm{ii}+\mathrm{jj}+\mathrm{kk})+$$
$$\frac{1}{2}(v_x w_x-v_y w_y)(\mathrm{ii}-\mathrm{jj})+$$
$$\frac{1}{6}(v_x w_x+v_y w_y-2v_z w_z)(\mathrm{ii}+\mathrm{jj}-2\mathrm{kk}).$$

To summarize, we've arrived at the following nine basis dyads:

$$\begin{array}{l} \mathsf{e}_1=(\mathrm{ii}+\mathrm{jj}+\mathrm{kk})/\sqrt{3} \\ \mathsf{e}_2=(\mathrm{ii}-\mathrm{jj})/\sqrt{2} \\ \mathsf{e}_3=(\mathrm{ii}+\mathrm{jj}-2\mathrm{kk})/\sqrt{6} \\ \mathsf{e}_4=(\mathrm{ij} +\mathrm{ji})/\sqrt{2} \\ \mathsf{e}_5=(\mathrm{jk} + \mathrm{kj})/\sqrt{2} \\ \mathsf{e}_6=(\mathrm{ki} +\mathrm{ik})/\sqrt{2} \\ \mathsf{e}_7=(\mathrm{ij} -\mathrm{ji})/\sqrt{2} \\ \mathsf{e}_8=(\mathrm{jk} - \mathrm{kj})/\sqrt{2} \\ \mathsf{e}_9=(\mathrm{ki} -\mathrm{ik})/\sqrt{2} \end{array} $$

In this basis, $\mathsf{D}=\mathrm{v}\mathrm{w}=D_i \mathsf{e}_i$ has the following components:

$$\begin{array}{l} D_1=(v_x w_x+v_y w_y+v_z w_z)/\sqrt{3} \\ D_2=(v_x w_x-v_y w_y)/\sqrt{2} \\ D_3=(v_x w_x+v_y w_y-2v_z w_z)/\sqrt{6} \\ D_4=(v_x w_y +v_y w_x)/\sqrt{2} \\ D_5=(v_y w_z + v_z w_y)/\sqrt{2} \\ D_6=(v_z w_x +v_x w_z)/\sqrt{2} \\ D_7=(v_x w_y -v_y w_x)/\sqrt{2} \\ D_8=(v_y w_z - v_z w_y)/\sqrt{2} \\ D_9=(v_z w_x -v_x w_z)/\sqrt{2}  \end{array} $$

As we will soon see, using this particular coordinate basis will greatly help us understand some of the fundamental, coordinate-independent properties of dyads.


### Operators

Unlike vectors, which have a straightforward physical interpretation as "arrows in space", dyads require some deeper context to understand what they are, and &mdash; perhaps more importantly &mdash; what they _do_.

Dyads can be regarded as _operators_. In fact, all tensors can be interpreted as operators. An operator is something that _does something_ to another tensor. For example, an operator acting on a physical vector can change that vector's magnitude, change its direction, or both.

Consider the case of changing a vector's magnitude. That can happen if we multiply the vector by a scalar:

$$\mathrm{v}\to m \mathrm{v}$$

In this example, the scalar ($m$) is acting as an operator. It's literally "scaling" the vector's magnitude by the factor $m$ &mdash; which explains the motivation behind the term "scalar" in the first place! Now, if $m$ has physical units, note that we necessarily end up with different physical units on the resulting vector compared to the one we started with. If $m$ is mass and $\mathrm{v}$ is a velocity vector, then the result $m\mathrm{v}$ is a momentum vector, usually labeled $\mathrm{p}$.

Vectors too can act as operators. The vector $\mathrm{v}$ can operate on the vector $\mathrm{w}$ to produce a scalar product, the value of which is defined as the product of the two vector magnitudes with an additional factor being the cosine of the angle $\alpha$ between the two vectors:

$$\mathrm{w}\to \mathrm{v}\cdot\mathrm{w}=|\mathrm{v}||\mathrm{w}| \cos \alpha$$ 

Note that the result of this "dot product" or "inner product" is coordinate independent.

Geometrically, the dot product between two vectors is one of vector's _projection_ onto the other vector multiplied by the length of that other vector. The result is the same regardless of which vector we project onto the other.

Instead of a physical vector, we can use a basis vector as the operator in the dot product. Due to the orthonormal nature of these basis vectors, they have the following dot products when operating amongst themselves:

$$\mathrm{i}\cdot\mathrm{i}=\mathrm{j}\cdot\mathrm{j}=\mathrm{k}\cdot\mathrm{k}=1$$
$$\mathrm{i}\cdot\mathrm{j}=\mathrm{j}\cdot\mathrm{i}=\mathrm{i}\cdot\mathrm{k}=\mathrm{k}\cdot\mathrm{i}=\mathrm{j}\cdot\mathrm{k}=\mathrm{k}\cdot\mathrm{j}=0$$

From the above, we can show that

$$\begin{array}{ccc}\mathrm{i}\cdot\mathrm{v}=v_x & \mathrm{j}\cdot\mathrm{v}=v_y & \mathrm{k}\cdot\mathrm{v}=v_z \end{array}$$

These results _are_ coordinate-dependent, which is not surprising given that the operators are themselves related to the coordinate basis. Relatedly, $\mathrm{v}\cdot\mathrm{w}$ can also be expressed as

$$\mathrm{v}\cdot\mathrm{w}=v_x w_x + v_y w_y + v_z w_z,$$

which, although still coordinate-independent, lacks the clear geometric interpretability that comes with the alternate expression $\mathrm{v}\cdot\mathrm{w}=|\mathrm{v}||\mathrm{w}| \cos \alpha$ we saw earlier.

The outer product we used to introduce dyads is another example of using a vector as an operator on another vector, in that case producing a dyad.

There are more ways vectors can behave as operators on other vectors, as we'll see later. 

### Dyads as operators

Earlier we remarked how, in general, a dyad can be representated as a _matrix_. In linear algebra, we learn that a fundamental use of matrices is to _transform_ vectors. The vector $\mathrm{x}$ can be transformed into the vector $\mathrm{y}$ using the matrix $\mathsf{A}$:

$$\mathrm{y}=\mathsf{A}\cdot \mathrm{x}$$

The details of how one calculates the product of $\mathsf{A}$ and $\mathrm{x}$ is generally presented as "matrix multiplication", an operation between arrays of coordinate-dependent components. The matrix would be a $3\times 3$ square array, while the vectors $\mathrm{x}$ and $\mathrm{y}$ would be "column" arrays of shape $3\times 1$.

If we instead represent $\mathsf{A}$ and $\mathrm{x}$ in terms of the basis vectors $\mathrm{i},\mathrm{j},\mathrm{k}$, and take the dot product operator to have the same meaning as before, the single dot product distributes out to a large number of dot products of basis elements, some of which survive and some of which vanish:

$$\mathsf{A}\cdot \mathrm{x}=$$
$$(A_1\mathrm{ii}+A_2\mathrm{ij}+...+A_9\mathrm{kk})\cdot(x_1\mathrm{i}+x_2\mathrm{j}+x_3\mathrm{k})=$$
$$A_1 x_1\mathrm{i}(\mathrm{i}\cdot\mathrm{i})+A_2 x_1\mathrm{i}(\mathrm{j}\cdot\mathrm{i})+...+A_9 x_3\mathrm{k}(\mathrm{k}\cdot\mathrm{k})=$$
$$(A_1 x_1+A_4 x_2+A_7 x_3)\mathrm{i}+$$
$$(A_2x_1+A_5 x_2+A_8 x_3)\mathrm{j}+$$
$$(A_3 x_1+A_6x_2+A_9 x_3)\mathrm{k}.$$

The result is a vector, equivalent to that given by "matrix multiplication", but here derived by applying the dot product distributed out onto all of the basis components.

But again, the above result applies to the highly specific context of using a particular coordinate basis.

### Geometric interpretation of dyad operators

We were able to interpret $\mathrm{v}\cdot\mathrm{w}$ in a purely geometric, coordinate-free way, so can we likewise develop a geometric, coordinate-free interpretation of $\mathsf{A}\cdot\mathrm{x}$?

In the previous section we saw how a dot product can "carve into" an outer product of vectors. Specifically, we can _associate_ the dot product operation in such a way that breaks the outer product up:

$$\mathsf{D}\cdot \mathrm{u}=\mathrm{vw}\cdot\mathrm{u}=\mathrm{v}(\mathrm{w}\cdot\mathrm{u})$$

In words, this is saying that the effect of applying the dyad $\mathsf{D}=\mathrm{vw}$ as an operator is a two-step process: first take the dot product with $\mathrm{w}$, then secondly scale $\mathrm{v}$ by the scalar result of the first step. This is a purely geometric, coordinate-free interpretation. 

The two-step process of taking a dot product and then scaling a vector is by no means the only operation we can apply to vectors using dyad operators. The caveat is that, in order to access these additional operations, we either need to sum together more than one outer-product dyad (such sums of dyads are referred to as _dyadics_), or, we can decompose a single dyad into components that embody these more pure operations. By "decompose into components" we specifically mean coordinate-free components, since we've already seen how to decompose a dyad into coordinate-dependent (basis) components like $\mathrm{ii}$, $\mathrm{ij}$, etc., as well as the $\mathsf{e}_i$ basis set we derived earlier.

In fact, the set of nine $\mathsf{e}_i$ basis elements from earlier is the key to framing the "pure" dyad operations. These pure operations are:

1. Isotropic scaling: the vector gets scaled by a value irrespective of its direction 

2. Asymmetric rotation: the vector changes its direction about an axis of rotation, with scaling in the plane perpendicular to that axis as it rotates 

3. Deviatoric scaling: relative to an orthonormal frame of three directions, vector components in up to two axes can scale independently and the third component scales in a compensating fashion


Let's now examine each of these in greater detail.

### Isotropic scaling

This operation is fully embodied by the first basis dyad (and its component) that we derived:

$$D_1\mathsf{e}_1=$$
$$\frac{1}{3}(v_x w_x+v_y w_y+v_z w_z)(\mathrm{ii}+\mathrm{jj}+\mathrm{kk})=\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})\mathsf{i}$$

There are two simplifications we've made here. 

First, the scalar component is literally just $1/3$ times the dot product $\mathrm{v}\cdot\mathrm{w}$. This shows that the dot (inner) product of two vectors is in fact just a subcomponent of the outer product of those vectors. 

The second simplification we made was in defining the _identity_ dyadic:

$$\mathsf{i}=\mathrm{ii}+\mathrm{jj}+\mathrm{kk}$$

which has the special property of mapping any vector onto itself:

$$\mathsf{i}\cdot\mathrm{v}=\mathrm{v}$$

The identity dyadic is often represented with an uppercase "I", but we've chosen lowercase here to help differentiate it with the moment of interia dyadic discussed later.

This is also a good place to introduce the shorthand notation $\mathrm{v}^2=\mathrm{vv}$ which denotes the outer product of a vector with itself. Using this notation the identity dyadic is written:

$$\mathsf{i}=\mathrm{i}^2+\mathrm{j}^2+\mathrm{k}^2$$

_Isotropic_ means "the same in all directions". This so-called _isotropic component_ of the dyad $\mathsf{D}=\mathrm{vw}$ represents an operator that (via dot producting) scales a vector $\mathrm{u}$ by the constant value $\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})$ irrespective of its direction:

$$\mathsf{D}_{\mathrm{iso}}\cdot\mathrm{u}=D_1\mathsf{e}_1\cdot\mathrm{u}=\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})\mathsf{i}\cdot\mathrm{u}=\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})\mathrm{u}$$

### Vector cross product

Besides the dot product and outer product, another way a vector $\mathrm{v}$ can operate on a vector $\mathrm{w}$ is the _cross product_. In Cartesian coordinates, the cross product is represented as:

$$\mathrm{v}\times\mathrm{w}=$$
$$(v_y w_z - v_z w_y)\mathrm{i}+(v_z w_x -v_x w_z)\mathrm{j}+(v_x w_y -v_y w_x)\mathrm{k}$$

which can be derived from the outer products of the Cartesian basis vectors:

$$\mathrm{i}\times\mathrm{i}=\mathrm{j}\times\mathrm{j}=\mathrm{k}\times\mathrm{k}=0$$
$$\mathrm{i}\times\mathrm{j}=\mathrm{k}\quad\mathrm{j}\times\mathrm{i}=-\mathrm{k}$$
$$\mathrm{j}\times\mathrm{k}=\mathrm{i}\quad\mathrm{k}\times\mathrm{j}=-\mathrm{i}$$
$$\mathrm{k}\times\mathrm{i}=\mathrm{j}\quad\mathrm{i}\times\mathrm{k}=-\mathrm{j}$$

Geometrically, the cross product $\mathrm{v}\times\mathrm{w}$ is a vector orthogonal to the plane spanned by $\mathrm{v}$ and $\mathrm{w}$, whose direction is determined by the "right hand rule" and whose magnitude is

$$|\mathrm{v}\times\mathrm{w}|=|\mathrm{v}||\mathrm{w}| \sin \alpha$$ 

where $\alpha$ is the angle between $\mathrm{v}$ and $\mathrm{w}$.

As it turns out, $|\mathrm{v}||\mathrm{w}| \sin \alpha$ is the _area_ of the 2-dimensional parallelogram bounded by the vectors $\mathrm{v}$ and $\mathrm{w}$. That being said, the vector product does _not_ contain all the detailed shape information of that parallelogram---that detailed shape information requires knowledge of _six_ elements, whereas the vector product has only three components (as do all vectors in three dimensions). Possible six-element sets include:

- the three components of each the two input vectors taken together,
- the cross-product's three components plus the three dot products $\mathrm{v}\cdot \mathrm{v}$, $\mathrm{v}\cdot\mathrm{w}$ and $\mathrm{w}\cdot\mathrm{w}$,

and there are more possible combinations.

In any event, the cross product can be interpreted geometrically as an oriented area element residing in the plane spanned by $\mathrm{v}$ and $\mathrm{w}$, the shape of that area not being encoded. In fact, _any_ vector can be interpreted in this manner, i.e., a directed area element instead of the usual directed line element.


### Asymmetric rotation

Now, the astute reader will have noticed that the components of $\mathrm{v}\times\mathrm{w}$ are the same as the $D_7,D_8,D_9$ dyad components we derived earlier. Specifically, one can see that:

$$D_7\mathsf{e}_7+D_8\mathsf{e}_8+D_9\mathsf{e}_9=\frac{1}{\sqrt{2}}(\mathsf{e}_7\mathrm{k}+\mathsf{e}_8\mathrm{i}+\mathsf{e}_9\mathrm{j})\cdot(\mathrm{v}\times\mathrm{w})$$

This quantity is referred to as the _asymmetric component_ of the dyad $\mathsf{D}=\mathrm{vw}$. Indeed, we had already associated $D_7,D_8,D_9$ with the asymmetric part of the dyad in our initial discussion of the $D_i\mathsf{e}_i$ basis decomposition.

The terms $\mathsf{e}_7\mathrm{k}$, $\mathsf{e}_8\mathrm{i}$, and $\mathsf{e}_9\mathrm{j}$ (and their sum) are our first encounter with _third-order tensors_, meaning tensors that have three indices. We also refer to these as _triadics_ (sums of triads), analogous to dyadics (sums of dyads). If we expand it out in terms of the basis vectors:

$$\frac{1}{\sqrt{2}}(\mathsf{e}_7\mathrm{k}+\mathsf{e}_8\mathrm{i}+\mathsf{e}_9\mathrm{j})=$$
$$\frac{1}{2}[(\mathrm{ij} -\mathrm{ji})\mathrm{k}+(\mathrm{jk} - \mathrm{kj})\mathrm{i}+(\mathrm{ki} -\mathrm{ik})\mathrm{j}]$$

The third-order tensor above in square brackets is known as the _Levi-Civita tensor_, represented as $\mathcal{E}$. We use caligraphic typeface for representing third-order tensors in an index-free manner.

To summarize:

$$D_7\mathsf{e}_7+D_8\mathsf{e}_8+D_9\mathsf{e}_9=\frac{1}{2} \mathcal{E}\cdot(\mathrm{v}\times \mathrm{w})$$

As an important aside: if we instead start with $\mathcal{E}$ and dot product it with two vectors $\mathrm{a}$ and $\mathrm{b}$ in sequence, we find:

$$(\mathcal{E}\cdot\mathrm{a})\cdot\mathrm{b}=$$
$$(a_y b_z - a_z b_y)\mathrm{i}+(a_z b_x -a_x b_z)\mathrm{j}+(a_x b_y -a_y b_x)\mathrm{k}$$
$$=\mathrm{a}\times\mathrm{b}$$

Taken together with the earlier equations, this reveals that dot-producting the antisymmetric component of a dyad $\mathsf{D}=\mathrm{vw}$ with a vector $\mathrm{u}$ is equivalent to:

$$\mathsf{D}_{\mathrm{anti}}\cdot\mathrm{u}=\frac{1}{2}(\mathcal{E}\cdot(\mathrm{v}\times\mathrm{w}))\cdot\mathrm{u}=\frac{1}{2}(\mathrm{v}\times\mathrm{w})\times\mathrm{u}$$

This operation _rotates_ the vector $\mathrm{u}$ by 90 degrees about the axis directed along $\mathrm{v}\times\mathrm{w}$, according to the right-hand rule, and scales it by factors $|\mathrm{v}\times\mathrm{w}|$ and the sine of the angle between $\mathrm{v}\times\mathrm{w}$ and $\mathrm{u}$. The scaling is isolated to the component of $\mathrm{u}$ which lies in the plane spanned by $\mathrm{v}$ and $\mathrm{w}$.

### Deviatoric scaling and eigenspaces

So far, we've decomposed the outer product $\mathrm{vw}$ into two geometrically meaningful components,

$$\mathsf{D}=\mathrm{vw}=\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})\mathsf{i} + \frac{1}{2}\mathcal{E}\cdot(\mathrm{v}\times \mathrm{w})+...$$
 
namely, an isotropic scaling component and an antisymmetric rotation component. 

The final component in our purely geometric decomposition, known as the _deviatoric scaling_ component, encapsulates the remaining five dyad basis elements $D_2$, $D_3$, $D_4$, $D_5$, and $D_6$.

Let's return to the $3\times 3$ matrix representation, noting that

$$
\begin{aligned}
\mathsf e_2 &= \frac{1}{\sqrt2}
\begin{pmatrix}
 1 & 0 & 0\\
 0 &-1 & 0\\
 0 & 0 & 0
\end{pmatrix}\quad
\mathsf e_3 = \frac{1}{\sqrt6}
\begin{pmatrix}
 1 & 0 & 0\\
 0 & 1 & 0\\
 0 & 0 &-2
\end{pmatrix}\\[6pt]
\mathsf e_4 &= \frac{1}{\sqrt2}
\begin{pmatrix}
 0 & 1 & 0\\
 1 & 0 & 0\\
 0 & 0 & 0
\end{pmatrix}\quad
\mathsf e_5 = \frac{1}{\sqrt2}
\begin{pmatrix}
 0 & 0 & 0\\
 0 & 0 & 1\\
 0 & 1 & 0
\end{pmatrix}\\[6pt]
\mathsf e_6 &= \frac{1}{\sqrt2}
\begin{pmatrix}
 0 & 0 & 1\\
 0 & 0 & 0\\
 1 & 0 & 0
\end{pmatrix}
\end{aligned}$$

From these, we can see that 

$$
\mathsf{D}_{\mathrm{dev}}=D_2\mathsf{e}_2+D_3\mathsf{e}_3+D_4\mathsf{e}_4+D_5\mathsf{e}_5+D_6\mathsf{e}_6
$$
$$
=\frac{1}{\sqrt2}
\begin{pmatrix}
 \dfrac{D_3}{\sqrt3}+ D_2 & D_4 & D_6\\[6pt]
 D_4 &  \dfrac{D_3}{\sqrt3}-D_2 & D_5\\[6pt]
 D_6 & D_5 & -\dfrac{2D_3}{\sqrt3}
\end{pmatrix}
$$

This is a _traceless symmetric_ matrix. Traceless means the diagonal components sum to zero. Even though it's represented as a 9-component matrix, there are still only 5 independent values that specify it.

But what does the deviatoric component actually do from a geometric operator standpoint? The key lies in the concept of an _eigenspace_. Specifically, any traceless symmetric $3\times3$ matrix can be decomposed into three orthonormal _eigenvectors_, and there are two independent _eigenvalues_.

Eigenvectors and eigenvalues carry a native geometric meaning in terms of their transformation of some vector $\mathrm{u}$. Each eigenvector specifies a direction along which the component of $\mathrm{u}$ is scaled by a constant value, namely the eigenvalue associated with that eigenvector. If the eigenvectors are denoted by $\mathrm{E}_1$, $\mathrm{E}_2$, $\mathrm{E}_3$ and the eigenvalues by $\lambda_1$, $\lambda_2$, $\lambda_3$, then this is expressed as:

$$\mathsf{D}_\mathrm{dev}\cdot\mathrm{u}=(\lambda_1\mathrm{E}_1\mathrm{E}_1+\lambda_2\mathrm{E}_2\mathrm{E}_2+\lambda_3\mathrm{E}_3\mathrm{E}_3)\cdot\mathrm{u}$$
$$=\lambda_1\mathrm{E}_1(\mathrm{E}_1\cdot\mathrm{u})+\lambda_2\mathrm{E}_2(\mathrm{E}_2\cdot\mathrm{u})+\lambda_3\mathrm{E}_3(\mathrm{E}_3\cdot\mathrm{u})$$

Note that the set of orthonormal eigenvectors completely span our three dimensions and therefore can construct the identity tensor:

$$\mathsf{i}=\mathrm{E}_1^2+\mathrm{E}_2^2+\mathrm{E}_3^2$$

The tracelessness property,

$$\lambda_1+\lambda_2+\lambda_3=0$$

is important because, geometrically, it means that the deviatoric operator can scale only up to two dimensions independently, the third dimension must be scaled in such a way to preserve the tracelessness. If we were to apply the deviatoric to a set of three or more vectors defining a _volume_, we'd find that the operator is _volume-preserving_. There is more to come on vector-derived volumes later in this text.

Unlike coordinate basis vectors, which are an arbitrarily defined external frame for our three dimensions, the eigenvectors are an orthonormal frame _intrinsic_ to the dyad $\mathsf{D}$. The eigenspace is just as intrinsic to the dyad as a magnitude is to a vector.

We _could_ technically compute the eigenspace for the $3\times 3$ matrix specified above, in terms of the $D_i$ components, using a particular algorithm for this purpose from linear algebra. However, that general closed form solution is very complicated and would not be easy to interpret or work with. Instead, we can reason what the eigenspace is from some key observations.

To begin, recall that the symmetric part of the dyad is

$$\mathsf{D}_\mathrm{symm}=\frac{1}{2}(\mathrm{vw}+\mathrm{wv})$$

Now, since the deviatoric is the symmetric part minus the $D_1\mathsf{e}_1$ (isotropic scaling) component,

$$\mathsf{D}_\mathrm{dev}=\frac{1}{2}(\mathrm{vw}+\mathrm{wv})-\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})\mathsf{i}$$

From here we can discover our first eigenvector, $\mathrm{E}_1$. _If_ $\mathrm{E}_1$ were orthogonal to both $\mathrm{v}$ and $\mathrm{w}$, _then_ $\mathrm{v}\cdot\mathrm{E}_1=\mathrm{w}\cdot\mathrm{E}_1=0$, and:

$$
\begin{aligned}
\mathsf{D}_\mathrm{dev}\cdot\mathrm{E}_1&=\frac{1}{2}(\mathrm{v}(\mathrm{w}\cdot\mathrm{E}_1)+\mathrm{w}(\mathrm{v}\cdot\mathrm{E}_1))-\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})(\mathsf{i}\cdot\mathrm{E}_1)\\
&=\frac{1}{2}(\mathrm{v}(0)+\mathrm{w}(0))-\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})(\mathrm{E}_1)\\
&=-\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})\mathrm{E}_1
\end{aligned}
$$

This says that $\mathsf{D}_\mathrm{dev}$ as an operator scales $\mathrm{E}_1$ by the constant $-\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})$. Therefore, $\mathrm{E}_1$ is _by definition_ an eigenvector of the deviatoric component of $\mathsf{D}$. Since this was based on the assumption that $\mathrm{E}_1$ is orthogonal to both $\mathrm{v}$ and $\mathrm{w}$, a natural choice (with normalization) is their cross product:

$$\mathrm{E}_1=\llbracket \mathrm{v}\times\mathrm{w} \rrbracket,$$

(though this doesn't work in the special case of $\mathrm v$ and $\mathrm w$ being collinear or antiparallel-- we'll address that later). We've also determined that the associated eigenvalue is

$$\lambda_1=-\frac{1}{3}(\mathrm{v}\cdot\mathrm{w}).$$

The other two eigenvectors must be orthogonal to $\mathrm{E}_1$ and therefore reside in the plane spanned by $\mathrm{v}$ and $\mathrm{w}$. This means they must be linear combinations of $\mathrm{v}$ and $\mathrm{w}$, i.e., 

$$
\mathrm{E}_2=a_2 \mathrm{v}+b_2 \mathrm{w} \\
\mathrm{E}_3=a_3 \mathrm{v}+b_3 \mathrm{w}
$$

and we want to solve for the coefficients $a_2,b_2,a_3,b_3$. Let's apply the deviatoric operator to a generic vector $\mathrm{x}=a \mathrm{v}+b\mathrm{w}$ in this plane space:

$$\mathsf{D}_\mathrm{dev}\cdot\mathrm{x}=\left(\frac{1}{2}(\mathrm{vw}+\mathrm{wv})-\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})\mathsf{i}\right)\cdot(a \mathrm{v}+b \mathrm{w})$$

Distributing out the dot product across all the terms eventually gives

$$\mathsf{D}_\mathrm{dev}\cdot\mathrm{x}=$$
$$\left(\frac{a (\mathrm{v}\cdot\mathrm{w})}{6}+\frac{b |\mathrm{w}|^2}{2}\right)\mathrm{v}+\left(\frac{b (\mathrm{v}\cdot\mathrm{w})}{6}+\frac{a |\mathrm{v}|^2}{2}\right)\mathrm{w}$$

For the special case where $\mathrm{x}$ is an eigenvector with eigenvalue $\lambda$, we know that 

$$\mathsf{D}_\mathrm{dev}\cdot\mathrm{x}=\lambda \mathrm{x}= \lambda (a \mathrm{v}+b \mathrm{w})$$

When combined with the previous vector equation, we obtain a system of two scalar equations (the components of $\mathrm{v}$ and $\mathrm{w}$):

$$
\lambda a =\frac{a (\mathrm{v}\cdot\mathrm{w})}{6}+\frac{b |\mathrm{w}|^2}{2}
$$
$$
\lambda b =\frac{b (\mathrm{v}\cdot\mathrm{w})}{6}+\frac{a |\mathrm{v}|^2}{2}
$$

Combining these by eliminating $\lambda$ and simplifying gives

$$
\left(\frac{a}{b}\right)^2=\left(\frac{|\mathrm{w}|}{|\mathrm{v}|}\right)^2
$$

Taking the square root of both sides,

$$\frac{a}{b}=\pm\frac{|\mathrm{w}|}{|\mathrm{v}|}$$

At this point we recognize the two possible values indicated by the "$\pm$" correspond to the two eigenvectors $\mathrm{E}_2$ and $\mathrm{E}_3$, i.e.,

$$\frac{a_2}{b_2}=\frac{|\mathrm{w}|}{|\mathrm{v}|} \quad\quad \frac{a_3}{b_3}=-\frac{|\mathrm{w}|}{|\mathrm{v}|}$$

Inserting these into the above equation with $\lambda$ yields the eigenvalues:

$$\lambda_2=\frac{\mathrm{v}\cdot\mathrm{w}}{6}+\frac{|\mathrm{v}||\mathrm{w}|}{2}\quad\quad \lambda_3=\frac{\mathrm{v}\cdot\mathrm{w}}{6}-\frac{|\mathrm{v}||\mathrm{w}|}{2}$$

Notice that $\lambda_1+\lambda_2+\lambda_3=0$, as expected since the deviatoric is _traceless_. 

Using the component-wise set of equations, we can solve to find

$$\begin{array}{cc} a_2=|\mathrm{w}| && b_2=|\mathrm{v}| \\ a_3=|\mathrm{w}| && b_3=-|\mathrm{v}| \end{array}$$

With normalization, the eigenvectors are therefore

$$\mathrm{E}_2=\llbracket |\mathrm{w}|\mathrm{v}+|\mathrm{v}|\mathrm{w}\rrbracket \quad\quad \mathrm{E}_3=\llbracket |\mathrm{w}|\mathrm{v}-|\mathrm{v}|\mathrm{w}\rrbracket $$

Bringing it all together, we can now fully expand the dyad into the isotropic scaling, antisymmetric rotational, and deviatoric symmetric components:

$$\begin{aligned}\mathsf{D}&=\mathrm{vw}\\ &=\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})\mathsf{i} + \frac{1}{2}\mathcal{E}\cdot(\mathrm{v}\times \mathrm{w})\\ &\quad -\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})\llbracket \mathrm{v}\times\mathrm{w}\rrbracket ^2 \\ &\quad+\left(\frac{\mathrm{v}\cdot\mathrm{w}}{6}+\frac{|\mathrm{v}||\mathrm{w}|}{2}\right)\llbracket |\mathrm{w}|\mathrm{v}+|\mathrm{v}|\mathrm{w}\rrbracket ^2 \\ &\quad+\left(\frac{\mathrm{v}\cdot\mathrm{w}}{6}-\frac{|\mathrm{v}||\mathrm{w}|}{2}\right)\llbracket |\mathrm{w}|\mathrm{v}-|\mathrm{v}|\mathrm{w}\rrbracket ^2\end{aligned}$$

### Special case of collinear vectors

Recall that when we were deriving the eigenvectors of the deviatoric symmetric part, we made an assumption that $\mathrm v$ and $\mathrm w$ have a defined cross-product, which is only true if they are neither collinear (pointing in the same direction) nor antiparallel (pointing in opposite directions). In the special cases of them being collinear or antiparallel, the deviatoric's eigenspace has a different structure and this will affect how we write the decomposition. 

If they are collinear or antiparallel, then we can actually find _two_ eigenvectors $\mathrm{E}'_1$ and $\mathrm{E}'_2$ which are orthogonal to both $\mathrm v$ and $\mathrm w$, using the same reasoning that led us to characterize $\mathrm{E}_1$ previously. As with our original $\mathrm{E}_1$, these two eigenvectors must have an eigenvalue of

$$\lambda'_1=\lambda'_2=-\frac{1}{3}(\mathrm{v}\cdot \mathrm{w})$$

The third eigenvector in this special case, $\mathrm{E}'_3$, is orthogonal to the first two eigenvectors and therefore collinear or antiparallel with $\mathrm v$ and $\mathrm w$. The normalization makes it straightforward to just use

$$\mathrm{E}'_3=\llbracket \mathrm{v} \rrbracket$$

and the eigenvalue is constrainted by the tracelessness property to be

$$\lambda'_3=\frac{2}{3}(\mathrm{v}\cdot\mathrm{w})$$

We don't actually have to define the $\mathrm{E}'_1$ and $\mathrm{E}'_2$ eigenvectors explicitly, because we can project onto the plane in which they reside via the operator

$$\mathsf{i}-\llbracket \mathrm{v} \rrbracket^2$$

and then scale by their common eigenvalue. Putting it all together, the deviatoric symmetric part for the case of collinear or antiparallel $\mathrm v$ and $\mathrm w$ is

$$\mathsf{D}_\mathrm{dev}=\frac{2}{3}(\mathrm{v}\cdot\mathrm{w})\llbracket \mathrm{v} \rrbracket^2-\frac{1}{3}(\mathrm{v}\cdot \mathrm{w})(\mathsf{i}-\llbracket \mathrm{v} \rrbracket^2)$$

However, notice that when we add in $\mathsf{D}_\mathrm{iso}=\frac{1}{3}(\mathrm{v}\cdot\mathrm{w})\mathsf{i}$ (also note the antisymmetric part is zero in this special case), we reduce to

$$\mathsf{D}=\mathrm{vw}=(\mathrm{v}\cdot\mathrm{w})\llbracket \mathrm{v}\rrbracket^2,$$

although keeping the isotropic and deviatoric parts separated is still instructive.

### Special case of orthogonal vectors

For the case where $\mathrm{v}$ and $\mathrm{w}$ are orthogonal, all terms that include $\mathrm{v}\cdot\mathrm{w}$ go to zero, and the expansion specializes to:

$$\mathsf{D}=\frac{1}{2}\mathcal{E}\cdot(\mathrm{v}\times \mathrm{w}) + \frac{|\mathrm{v}||\mathrm{w}|}{2}\left(\llbracket |\mathrm{w}|\mathrm{v}+|\mathrm{v}|\mathrm{w}\rrbracket^2 - \llbracket |\mathrm{w}|\mathrm{v}-|\mathrm{v}|\mathrm{w}\rrbracket^2\right).$$

Note that the isotropic term in particular has disappeared entrely. Furthermore, since for the orthogonal case we have $|\mathrm{v}\times \mathrm{w}|=|\mathrm{v}||\mathrm{w}|$, we can pull out the scaling factor:

$$\mathsf{D}=\frac{|\mathrm{v}||\mathrm{w}|}{2}\left(\mathcal{E}\cdot\llbracket\mathrm{v}\times \mathrm{w}\rrbracket + \llbracket |\mathrm{w}|\mathrm{v}+|\mathrm{v}|\mathrm{w}\rrbracket^2 - \llbracket |\mathrm{w}|\mathrm{v}-|\mathrm{v}|\mathrm{w}\rrbracket^2\right).$$

Collectively the three terms serve to map (project and scale) any vector that this operates on into the plane spanned by $\mathrm{v}$ and $\mathrm{w}$. Indeed, as is always true for the dyad $\mathrm{vw}$, the final result is always parallel to v, since

$$\mathsf{D}\cdot \mathrm{u}=\mathrm{vw}\cdot\mathrm{u}=\mathrm{v}(\mathrm{w}\cdot\mathrm{u}),$$

even though the above decomposition obscures this.

Moreover, the following three unit vectors form an orthonormal basis:

$$\mathrm{e}_1=\llbracket\mathrm{v}\times \mathrm{w}\rrbracket$$
$$\mathrm{e}_2=\llbracket |\mathrm{w}|\mathrm{v}+|\mathrm{v}|\mathrm{w}\rrbracket$$
$$\mathrm{e}_3=\llbracket |\mathrm{w}|\mathrm{v}-|\mathrm{v}|\mathrm{w}\rrbracket,$$

in terms of which the decomposition becomes the rather compact:

$$\mathsf{D}=\frac{|\mathrm{v}||\mathrm{w}|}{2}\left(\mathcal{E}\cdot\mathrm{e}_1 + \mathrm{e}_2^2-\mathrm{e}_3^2\right).$$

### Matrix multiplication

So far we've taken the approach of trying to understand the dyad $\mathsf{D}=\mathrm{vw}$ in terms of geometrically meaningful components and with regard to how it behaves as an operator. A different perspective is to instead examine what happens when we operate on the dyad itself with a different operator.

Consider the operator $\mathsf{D}=\mathrm{vw}$ transforming the vector $\mathrm{u}$ into the vector $\mathrm{y}=\mathsf{D}\cdot\mathrm{u}$. If we simply append (via outer product) another vector $\mathrm{x}$ to the right side, this additional vector doesn't participate in the dot product contraction, only the vectors directly adjacent to the dot operator are affected:

$$(\mathsf{D}\cdot\mathrm{u})\mathrm{x}=\mathrm{vw}\cdot\mathrm{ux}=\mathrm{v}(\mathrm{w}\cdot\mathrm{u})\mathrm{x}=(\mathrm{v}\mathrm{w}\cdot\mathrm{u})\mathrm{x}=\mathrm{yx}.$$

Defining $\mathsf{F}=\mathrm{ux}$, and $\mathsf{G}=\mathrm{yx}$, this becomes

$$\mathsf{D}\cdot\mathsf{F}=\mathsf{G}.$$

This is another example of the "matrix multiplication" concept we discussed earlier, except this time between two dyads. In Einstein notation:

$$D_{ij}F_{jk}=G_{ik}$$

The key takeaway is that the dot operator contracts the _left tensor's right index with the right tensor's left index_. Simply writing $\mathsf{D}\cdot\mathsf{F}$ doesn't make this obvious. For instance, one might incorrectly assume they contract on both left ($D_{ji}F_{jk}$), both right ($D_{ij}F_{kj}$), or on their outside indices ($D_{ji}F_{kj}$). Those other contractions produce different results in general, hence the meaning of the dot product is important.

### Double dot product

Just as the regular dot product denotes contraction of a single index between two tensors, the double dot product "$:$" denotes a contraction between two indices:

$$\mathsf{A}:\mathsf{B}=A_{ij}B_{ji}$$

In exactly the same sense that the implied index contraction was important to acknowledge for $\mathsf{D}\cdot\mathsf{F}$, here we again require a specific contraction order. The rule is to contract the furthest right index belonging to the tensor on the left of the dot operation with the furthest left index belonging to the tensor on the right of the dot operation, then continue working outward until all the "dots" have been accounted for. We'll refer to this as "inside-outward" contraction, and it has a characteristic reflected ordering of the contracted indices on either side of the contraction operator (observe this in $A_{ij}B_{ji}$).

Taking the identity tensor $\mathsf{i}$ and the Levi-Civita tensor $\mathcal{E}$ as operators, the double dot product provides some important relationships:

$$\mathsf{i}:\mathrm{vw}=\mathrm{v}\cdot\mathrm{w}$$
$$\mathcal{E}:\mathrm{vw}=\mathrm{v}\times\mathrm{w}$$

We already had discovered previously that the dot product $\mathrm{v}\cdot\mathrm{w}$ and the cross product $\mathrm{v}\times\mathrm{w}$ are intrinsic to the isotropic and antisymmetric components of the dyad, so it shouldn't be too surprising that we can pull these entities out. It is striking however that they emerge with such simple operations.

A few more interesting relationships produced by the double dot product are:

$$\mathcal{E}:\mathsf{i}=\mathsf{i}:\mathcal{E}=0$$

$$\mathcal{E}:\mathcal{E}=-2\mathsf{i}$$

Lastly, the _trace_ of a matrix $\mathsf{A}$, a notable subject in linear algebra, is also expressable with the double dot product:

$$\mathrm{trace}(\mathsf{A})=\mathsf{i}:\mathsf{A}.$$

Recall the concept of trace arose earlier in the context of the _traceless_ symmetric deviatoric operator.

### Antisymmetrizer and related four-index tensors

Earlier we saw that the antisymmetric part of the dyad $\mathsf{D}=\mathrm{vw}$ can be written

$$\mathsf{D}_\mathrm{anti}=\frac{1}{2}(\mathrm{vw}-\mathrm{wv})=\frac{1}{2}\mathcal{E}\cdot(\mathrm{v}\times\mathrm{w}).$$

By substituting the equivalent expression involving the double dot product, this takes the form:

$$\mathsf{D}_\mathrm{anti}=\frac{1}{2}\mathcal{E}\cdot(\mathcal{E}:\mathrm{v}\mathrm{w})=\left(\frac{\mathcal{E}\cdot\mathcal{E}}{2}\right):\mathrm{v}\mathrm{w}=\mathfrak{A}:\mathrm{v}\mathrm{w},$$

where we have defined a four-index tensor

$$\mathfrak{A}=\frac{1}{2}\mathcal{E}\cdot\mathcal{E}.$$

This is known as the _antisymmetrizer_ tensor. We use Fraktur typeface to represent four-index tensors in the index-free manner. Just to make the implied contraction policy extra clear, in Einstein notation this is:

$$\mathfrak{A}\_{ijkl}=\frac{1}{2}\mathcal{E}\_{ijm}\mathcal{E}\_{mkl}.$$

The antisymmetrizer tensor gives the antisymmetric part of _any_ general dyadic $\mathsf{A}$, not just pure dyads:

$$\mathsf{A}_\mathrm{anti}=\mathfrak{A}:\mathsf{A}.$$

There are a few other four-index tensors that relate to the antisymmetrizer. The first is the _symmetrizer_ $\mathfrak{S}$, whose function can be directly inferred based on the earlier statement that any dyadic is the sum of its symmetric and antisymmetric parts:

$$\mathsf{A}_\mathrm{symm}=\mathfrak{S}:\mathsf{A}=\mathsf{A}-\mathfrak{A}:\mathsf{A}.$$

We can "factor out" the $\mathsf{A}$ in the right side of the above equation by introducing a four-index identity tensor $\mathfrak{I}$ which satisfies $\mathfrak{I}:\mathsf{A}=\mathsf{A}$ for any $\mathsf{A}$:

$$\mathfrak{S}:\mathsf{A}=\mathsf{A}-\mathfrak{A}:\mathsf{A}=(\mathfrak{I}-\mathfrak{A}):\mathsf{A},$$

or more succinctly:

$$\mathfrak{S}=\mathfrak{I}-\mathfrak{A}.$$

Note that $\mathfrak{I}$ is an identity operator only when used with the double dot product, analogous to how $\mathsf{i}$ is an identity operator only when used with the single dot product.

The last four-index tensor in this set pertains to the fact that the symmetric part minus the antisymmetric part is equal to the transpose:

$$\mathsf{A}^{\mathsf{T}}=\mathfrak{S}:\mathsf{A}-\mathfrak{A}:\mathsf{A}=(\mathfrak{S}-\mathfrak{A}):\mathsf{A}=\mathfrak{T}:\mathsf{A},$$

or equivalently,

$$\mathfrak{T}=\mathfrak{S}-\mathfrak{A},$$

which is the _transposer_ tensor.

Among these four-index tensors, we've thus far only expressed $\mathfrak{A}$ in terms of lower-index tensors (namely, in terms of the three-index $\mathcal{E}$). But all of them can be written in terms of unit-basis-vector outer products in the following manner:

$$\mathfrak{A}=\frac{1}{2}[\mathrm{ij}(\mathrm{ji}-\mathrm{ij})+\mathrm{jk}(\mathrm{kj}-\mathrm{jk})+\mathrm{ki}(\mathrm{ik}-\mathrm{ki})+\cdots]$$

$$\mathfrak{S}=\mathrm{iiii}+\mathrm{jjjj}+\mathrm{kkkk}+\frac{1}{2}[\mathrm{ij}(\mathrm{ji}+\mathrm{ij})+\cdots)]$$

$$\mathfrak{I}=\mathrm{iiii}+\mathrm{ijji}+\mathrm{ikki}+\cdots$$

$$\mathfrak{T}=\mathrm{iiii}+\mathrm{ijij}+\mathrm{ikik}+\cdots$$

One can show that $\mathfrak{A}$ has 12 terms, $\mathfrak{S}$ has 15 terms, $\mathfrak{I}$ and $\mathfrak{T}$ have 9 terms.

### Dyad operations revisited with four-index tensors

Previously we showed that a dyad expands out into three distinct, pure operations with physical interpretation: isotropic scaling, asymmetric rotation, and deviatoric scaling. The double dot product along with four-index tensors are a convenient way of specializing or filtering the original dyad operator into one or more of these pure operations. We already noted this for the antisymmetric operation:

$$\mathsf{D}_\mathrm{anti}=\mathfrak{A}:\mathsf{D}=\mathfrak{A}:\mathrm{vw}=\frac{1}{2}\mathcal{E}\cdot(\mathrm{v}\times\mathrm{w}).$$

Expanding upon the earlier interpretation that $\mathsf{D}$ is an operator which can operate on the vector $\mathrm{u}$, this operation being written as $\mathsf{D}\cdot\mathrm{u}$, then we can similarly write

$$(\mathfrak{A}:\mathsf{D})\cdot\mathrm{u}=\frac{1}{2}(\mathrm{v}\times\mathrm{w})\times\mathrm{u}$$

where $(\mathfrak{A}:\mathsf{D})$ is the effective operator formed by altering the original operator $\mathsf{D}$.

For the symmetric part, we can separate out the isotropic (non-zero trace) part from the deviatoric (traceless) part:

$$\mathsf{D}_{\mathrm{symm}}=\mathfrak{S}:\mathsf{D}=\frac{1}{3}\mathsf{i}\mathsf{i}:\mathsf{D}+\left(\mathfrak{S}-\frac{1}{3}\mathsf{i}\mathsf{i}\right):\mathsf{D}$$

$$\mathsf{D}_{\mathrm{iso}}=\left(\frac{\mathsf{i}^2}{3}\right):\mathsf{D}$$

$$\mathsf{D}_{\mathrm{dev}}=\left(\mathfrak{S}-\frac{\mathsf{i}^2}{3}\right):\mathsf{D}=\mathfrak{D}:\mathsf{D}$$

Here we used the above-noted $\mathrm{trace}(\mathsf{D})=\mathsf{i}:\mathsf{D}$, and for the specific case of $\mathsf{D}=\mathrm{vw}$, we have $\mathrm{trace}(\mathsf{D})=\mathrm{v}\cdot\mathrm{w}$. We also defined another four-index tensor

$$\mathfrak{D}=\mathfrak{S}-\frac{\mathsf{i}^2}{3}$$

which produces the deviatoric part of the dyadic it operates on (and happens to have 21 terms). Note that $\mathsf{i}^2$ is another four-index tensor (having 9 terms), we simply aren't giving it a dedicated label like the others.


### Four-dot product and trace of the matrix product

The outer product of two dyadics, $\mathsf{M}$ and $\mathsf{N}$, is written $\mathsf{MN}$. This outer product is a four-index tensor. Defining the four-dot product "$:\hspace{0.2em}:$" as the inside-outward contraction on four indices, there is an interesting relationship between the four-tensor $\mathfrak{I}$ and the matrix multiplication product ($\mathsf{M}\cdot\mathsf{N}$): 

$$\mathfrak{I}:\hspace{0.2em}:\mathsf{MN}=\mathsf{i}:(\mathsf{M}\cdot\mathsf{N})=\mathrm{trace}(\mathsf{M}\cdot\mathsf{N})=M\_{ij}N\_{ji}.$$

The fact that the unit basis terms in $\mathfrak{I}$ have the same $abba$ pattern seen in the contracted indices of $M_{ij}N_{ji}$ is an intuitive confirmation of this relationship. 

Note this is similar to the earlier-observed relationship between $\mathsf{i}$ and the vector inner product $\mathsf{v}\cdot\mathsf{w}$:

$$\mathsf{i}:\mathrm{vw}=\mathrm{v}\cdot\mathrm{w}$$

### Determinant: geometric concept

Consider the triad $\mathcal{V}$ composed of the triple outer product of the three vectors $\mathrm{v}$, $\mathrm{w}$, and $\mathrm{u}$:

$$\mathcal{V}=\mathrm{vwu}$$

Provided none of the three vectors are collinear, they span a volume $V$ given by:

$$V=\mathcal{E}\hspace{0.2em}\vdots\hspace{0.2em}\mathrm{vwu}=(\mathrm{v}\times\mathrm{w})\cdot\mathrm{u}$$

where we have introducted the triple dot product "$\hspace{0.2em}\vdots\hspace{0.2em}$", which conveys the same inside-outward index contraction policy that we saw for the double and single dot products. In Einstein notation:

$$V=\mathcal{E}\_{ijk}v\_k w\_j u\_i$$

which shows the characteristic reflected ordering of contracted indices on either side of the operator.

For an orthonormal basis, such as the usual $\mathrm{i,j,k}$, the volume is unity:

$$\mathcal{E}\hspace{0.2em}\vdots\hspace{0.2em}\mathrm{ijk}=1$$

Now if we take some general dyadic operator $\mathsf{M}$ and apply it to all the unit vectors, such that they are transformed and form a new transformed volume, the measure of that transformed volume is defined as the _determinant_ of $\mathsf{M}$:

$$\mathrm{det}(\mathsf{M})=\mathcal{E}\hspace{0.2em}\vdots\hspace{0.2em}(\mathsf{M}\cdot\mathrm{i})(\mathsf{M}\cdot\mathrm{j})(\mathsf{M}\cdot\mathrm{k})$$

Though well-defined, this is _not_ a coordinate-independent definition for the determinant, because it requires explicitly including the $\mathrm{i,j,k}$ basis vectors. In the same sense that we can write $\mathrm{trace}(\mathsf{M})=\mathsf{i}:\mathsf{M}$, we'd like to be able to write

$$\mathrm{det}(\mathsf{M})=\text{Д}\hspace{0.2em}\vdots\hspace{0.2em}\vdots\hspace{0.2em}\mathsf{M}^3$$

where Д is some yet-to-be-constructed six-index tensor (and the Cyrillic character for "D", as we will use Cyrillic characters to represent six-index tensors). Here we also have introduced the six-dot product "$\hspace{0.2em}\vdots\hspace{0.2em}\vdots\hspace{0.2em}$", which just as the other $n$-dot products, follows the inside-outward contraction, represented in Einstein notation as

$$\text{Д}\hspace{0.2em}\vdots\hspace{0.2em}\vdots\hspace{0.2em}\mathsf{M}^3=\text{Д}\_{ijklmn}M\_{nm}M\_{lk}M\_{ji}$$

again having the reflected ordering of contracted indices on either side of the dot operator. 

### Matrix multiplication from a six tensor

Whenever a operator tensor operates on an input tensor to produce an output tensor, the following relationship must hold: The number of indices on the output tensor is equal to the number of indices on the input tensor plus the number of indices on the operator tensor, minus twice the number of contraction dots involved in the operation.

From this, we can surmise that there ought to exist a six tensor, Ж, which, when operating via four-dot product on the outer product of two dyadics, $\mathsf{MN}$, produces their two-index matrix product, $\mathsf{M}\cdot \mathsf{N}$:

$$Ж:\hspace{0.2em}:\mathsf{M}\mathsf{N}=\mathsf{M}\cdot\mathsf{N}.$$

This somewhat resembles a relation we saw earlier,

$$\mathfrak{I}:\hspace{0.2em}:\mathsf{MN}=\mathrm{trace}(\mathsf{M}\cdot\mathsf{N}),$$

except for the fact that now we have two free (uncontracted) indices on the output instead of them all having been contracted.

Thinking in terms of the basis element patterns in these inside-outward contractions: Uncontracted indices (and corresponding constituent basis vectors) on the left side of the operator can correspond to contracted indices on the right side of the operator. Consider the following cases we've already seen:

$$\mathfrak{I}:\mathsf{M}=\underbrace{(ab}\_{\text{free}}\underbrace{ba)}\_{\text{contracted}}:\mathsf{M}=\mathsf{M}$$

$$\mathfrak{T}:\mathsf{M}=\underbrace{(ab}\_{\text{free}}\underbrace{ab)}\_{\text{contracted}}:\mathsf{M}=\mathsf{M}^\mathsf{T}$$

We can reason that this latest case adds a third freely-varying index, $c$, but only in the contracted portion of the operator, and in such a way that produces the _nested_ contraction in the required position:

$$Ж:\hspace{0.2em}:\mathsf{MN}=\underbrace{(ab}\_{\text{free}}\underbrace{bcca)}\_{\text{contracted}}:\hspace{0.2em}:\mathsf{MN}=\mathsf{M}\cdot\mathsf{N}$$

Just as $abba$ and $abab$ are the patterns for $\mathfrak{I}$ and $\mathfrak{T}$, the pattern for Ж is $abbcca$, meaning it has terms like:

$$Ж=\mathrm{ijjkki}+\mathrm{jkkiij}+\mathrm{kiijjk}+\cdots$$

Note that the example terms shown here exemplify the $abbcca$ pattern in cases where $a\ne b\ne c$. But Ж also contains terms for which $a=b\ne c$, for which $a\ne b= c$, for which $a=c\ne b$, and for which $a=b=c$. So we can equally write

$$Ж=\mathrm{ijjiii}+\mathrm{jjjkkj}+\mathrm{iiiiii}+\cdots$$

and the full expression for Ж has 27 total terms.

### Construction of the determinant operator

For a matrix with elements

$$\mathsf{A}=\begin{pmatrix}
 a & d & g\\
 b & e & h\\
 c & f & i
\end{pmatrix}$$

the formula for the determinant (given by, _e.g._, the "rule of Sarrus" in linear algebra) is

$$\text{det}(\mathsf{A})=aei+dhc+gbf-ceg-fha-ibd$$

Mapping these terms and their factors into the dyad basis elements, we can write out a possible form for Д directly in such a way that is consistent with the inside-outward contraction. The dyad basis element for $a$ is $\mathrm{ii}$, $b$ is $\mathrm{ji}$, and so on, yielding in its entirety:

$$\text{Д}\stackrel{?}{=}\mathrm{iijjkk}+\mathrm{jikjik}+\mathrm{kiijjk}-\mathrm{ikjjki}-\mathrm{jkkjii}-\mathrm{kkijji}$$

The reason we used "$\stackrel{?}{=}$" instead of "$=$" is that, within each term, we are free to rearrange the dyad basis elements, because the resulting six-dot product of that term's three dyad factors with $\mathsf{M}^3$ will be the same irrespective of their ordering. These reorderings however are not strictly equal to one another prior to any contraction, so they cannot all equal Д.

Examining the terms in our above proposed form for Д, we identify some key patterns:

1. The first term, $\mathrm{iijjkk}$, is one of the 27 terms that appear in $\mathsf{i}^3$.

2. The second term, $\mathrm{jikjik}$, is characterized by a _cyclic_ pattern $abcabc$ wherein the first three factors are verbatim repeated as the second set of three factors. Define the six-index tensor Ш as the 27-term sum of all such possible $abcabc$ combinations:

$$\text{Ш}=\mathrm{ijkijk}+\mathrm{jikjik}+\mathrm{kijkij}+\mathrm{kjikji}+\mathrm{ikjikj}+\mathrm{jkijki}+\cdots$$

3. The third term, $\mathrm{kiijjk}$, is one of the 27 terms found in the six tensor Ж we explored in the previous section which gives matrix multiplication.

4. The fourth term, $\mathrm{ikjjki}$, is one belonging to the "reflected" $abccba$ pattern, making it a member of the six-tensor/triple-dot-product equivalent of the identity tensor (analogous to $\mathfrak{I}$ for double dot product and $\mathsf{i}$ for single dot product). Call this six-tensor identity Б, which for any three-index tensor $\mathcal{F}$ has the behavior

$$Б\hspace{0.2em}\vdots\hspace{0.2em}\mathcal{F}=\mathcal{F}$$

<p style="margin-left: 2.5em">Just as Б is analogous to $\mathfrak{I}$, note that the Ш we defined earlier is analogous to $\mathfrak{T}$. Whereas $\mathfrak{T}$ swaps the two indices of a dyadic, Ш swaps the first and third indices on a triadic (without changing the middle second index). Also note that both of these have a similar sign-changing effect on $\mathcal{E}$:

$$\mathfrak{T}:\mathcal{E}=-\mathcal{E}\qquad Ш\hspace{0.2em}\vdots\hspace{0.2em}\mathcal{E}=-\mathcal{E}$$
</p>

5. The fifth term is a member of $\mathfrak{I}\mathsf{i}$.

6. The sixth term is a member of $\mathsf{i}\mathfrak{I}$.

Given the above observations, we might expect to be able to write Д as some linear combination of $\mathsf{i}^3$, Ш, Ж, Б, $\mathfrak{I}\mathsf{i}$, and $\mathsf{i}\mathfrak{I}$.

Though an exhaustive enumeration and comparison of the terms, many of which cancel out or stack up as multiples, we ultimately find:

$$Д=\frac{1}{6}(\mathsf{i}^3+Ш+Ж-Б-\mathfrak{I}\mathsf{i}-\mathsf{i}\mathfrak{I})$$

The earlier form of Д we proposed above, based on the Sarrus rule, is a special case of this result obtained by combining terms that are equivalent through the dyad reordering.

### Contributing terms in the determinant

Now that we have a suitable form for Д, we can evaluate the terms contributing to the determinant in their coordinate-free index-free form. Distribuing out the six-dot operator:

$$\begin{aligned}\text{det}(\mathsf{M})=Д\hspace{0.2em}\vdots\hspace{0.2em}\vdots\hspace{0.2em}\mathsf{M}^3=\frac{1}{6}(&\mathsf{i}^3\hspace{0.2em}\vdots\hspace{0.2em}\vdots\hspace{0.2em}\mathsf{M}^3+Ш\hspace{0.2em}\vdots\hspace{0.2em}\vdots\hspace{0.2em}\mathsf{M}^3+Ж\hspace{0.2em}\vdots\hspace{0.2em}\vdots\hspace{0.2em}\mathsf{M}^3 \\
&-Б\hspace{0.2em}\vdots\hspace{0.2em}\vdots\hspace{0.2em}\mathsf{M}^3-\mathfrak{I}\mathsf{i}\hspace{0.2em}\vdots\hspace{0.2em}\vdots\hspace{0.2em}\mathsf{M}^3-\mathsf{i}\mathfrak{I\hspace{0.2em}\vdots\hspace{0.2em}\vdots\hspace{0.2em}\mathsf{M}^3})\end{aligned}$$

Let's look at each term:

1. $\mathsf{i}^3\hspace{0.2em}\vdots\hspace{0.2em}\vdots\hspace{0.2em}\mathsf{M}^3$ --- Each of the three $\mathsf{i}$ instances cleanly annihilates with one of the three $\mathsf{M}$ instances to produce a factor $\mathsf{i}:\mathsf{M}=\text{trace}(\mathsf{M})$. So this term overall is $\text{trace}(\mathsf{M})^3$.