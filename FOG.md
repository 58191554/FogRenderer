# Fog Renderer

Zhen Tong 2023/7/17

Based on the ray tracer, the fog renderer is implemented by two part, the fog particle intersection and fog scattering. 

## Particle Intersection

First, we assume the fog particle is uniformly distributed in the atmosphere. The ray will hit fog particles before intersect with the object. To define when the ray hit the particle, exponential distribution is used to describe the probability of time interval between two event (ray hit particle). Taking the inverse of exponential distribution CDF, we can sample the fog with the same sample rate of no fog, and reach a fog-blur effect. 

<img src="C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230717121626474.png" alt="image-20230717121626474" style="zoom:50%;" />

According to the definition of the Poisson process, the time interval between the $k^{th}$ random event and the $k+1^{th}$ random event follows an exponential distribution [Eric P. Lafortune and Yves D. Willems]. And the probability of zero random event occurring in a time period of length $t$ is equal to
$$
P(t|hit=0)=\frac{e^{-\lambda t}(\lambda t)^0}{0!} = e^{-\lambda t}\\
P(t|hit>0) = 1-e^{-\lambda t}\\
\text{CDF of exponential distribution:}
F(t) = \begin{cases}
1 - e^{-\lambda t} & \text{if } t \geq 0 \\
0 & \text{if } t < 0
\end{cases}
$$
Then use uniform random variable $z\sim U[0, 1]$, we can sample particle intersection $t = -\frac{\log{(1-z)}}{\lambda}$, the physical interpretation of $\lambda$ is number of events per unit time, therefore, larger $\lambda$ will have shorter time $t$. 

| $\lambda=0.002$                                              | $\lambda=0.02$                                               | $0.2$                                                        |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| ![image-20230717130335139](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230717130335139.png) | ![image-20230717131028918](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230717131028918.png) | ![image-20230717132254982](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230717132254982.png) |

There is an assumption in the ray transmitting in the fag media, the radiance will go through absorbing, emitting, out-scattering, and in-scattering on its way. When we simplify the model to only focus on the absorbing and out-scattering, two process to decay the radiance in the model. According to Wojciech Jarosz: 
$$
dL = -\sigma_t(t) Ldt\\
L = L_0 e^{-\int_0^{t_{end}}\sigma_{t_{}}(x)dx} = L_0e^{-\sigma_tt_{end}}
$$
$dL$ is the light intensity derivative, and $dt$ is the time derivative along the ray. $\sigma_t=\sigma_a+\sigma_s$ is the sum of the absorption coefficient and the out-scattering coefficient. The first equation describes the decay of a the intensity $L$ when the ray go through place $\vec{o}+\vec{d}t$ the absorbing and out-scattering media. Solving the equation, we get the integrated formula that can tell us how much light after decay we observe by the intersection $t_{end}$.

## Fog Scattering

After the fog particle intersection, we need to implement fog particle reflect the light. By Monte Carlo sampling, we can have different probability among directions. According to the Henyey-Greenstein phase function,
$$
P_{HG}(\theta) = \frac{1}{4\pi}\frac{1-g^2}{(1+g^2-2g\cos{\theta})^{\frac{3}{2}}}
$$
Setting different value of g, we can probability distribution with $\theta$, the ray direction is from the negative side of polar axis. When g = 0, the BSDF is diffuse like; when g is negative, the ray transmit forward; and when g is positive, the ray reflect backward.

 ![image-20230717200549794](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230717200549794.png)

With different g, we can render different effect of fog. From the show case below, we can see the fog particle below the light is observed, while the fog particle far from the light is less visible since there is few reflection to the camera.

| g = 0                                                        | g=-0.7                                                       | g=0.8                                                        |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| ![image-20230717145248014](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230717145248014.png) | <img src="C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230717144404549.png" alt="image-20230717144404549" style="zoom:100%;" /> | <img src="C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230717144059903.png" alt="image-20230717144059903" style="zoom:100%;" /> |



## Challenge

1. The fog reflection by Henyey-Greenstein phase function has some problem, there isn't much different between g=-0.7 and g=0.8, and when g=0, it is more like the situation when g=0.8.

2. To better illustrate the fog effect, I want to use the spot light instead of the area light or point light. But, the effect of spot light cone is not good, though I adjust the g value to -1, this is probably because the reflection process is incorrect. You can see there is no light cone in my renderer.

   | My Renderer                                                  | Blender                                                      |
   | ------------------------------------------------------------ | ------------------------------------------------------------ |
   | ![image-20230717200744365](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230717200744365.png) | <img src="C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230717205856987.png" alt="image-20230717205856987" style="zoom:50%;" /> |

   

## Next Step

1. Solve the problem issued above
2. Implement mirror and glass material BSDF to build a better render.

## Reference

1. Light Transport in Participating Media: https://cs.dartmouth.edu/~wjarosz/publications/dissertation/chapter4.pdf
2. Rendering Participating Media with Bidirectional Path Tracing: http://luthuli.cs.uiuc.edu/~daf/courses/rendering/papers/lafortune96rendering.pdf
3. On Sampling Of Scattering Phase Functions: https://arxiv.org/pdf/1812.00799.pdf
4. Volumetric Path Tracing: https://justgood.dev/docs/740-paper.pdf

5. Blogs: 

   https://zhuanlan.zhihu.com/p/526598533

   https://zhuanlan.zhihu.com/p/56710440







![image-20230727113820683](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230727113820683.png)

![image-20230727115857241](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230727115857241.png)

spot light can not use random sample procedure, because the light cannot precisely point to the spot point position. 

​	

![image-20230727173222074](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230727173222074.png)

![image-20230727173249423](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230727173249423.png)