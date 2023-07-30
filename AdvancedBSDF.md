# Advanced BSDF Zhen TONG

## Mirror and Glass Material

### Mirror Material

Mirror implementation is straightforward. When the ray from the camera hit the mirror-surface object, the reflected ray is determined without Monte Carlo sample. The only thing need to do is emit a new ray from the hit point, and render the image with at least 2 ray depth, because the mirror reflection uses one depth, and the next object the reflected ray hits takes another ray depth. We don't need to multiply the $cos(\theta)$ here, because mirror surface does not follow the Lambertian Law.

| Ray Depth | 1                                                            | 2                                                            | 3                                                            |
| --------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| **Image** | ![image-20230724115041943](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230724115041943.png) | ![image-20230724115049584](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230724115049584.png) | ![image-20230724120612249](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230724120612249.png) |

As the picture shown above, the more ray depth we set, the rays on the dragon's body that did not hit the surface of the radiance (black) after the first reflection hit the surface of the radiance during the second round of reflection.

### Refraction

When the ray shot into some transparent media, like water and glass, the ray will refract instead of reflect. In this case, we need apply the Snell's rule. To represent the direction vector of a ray with two angle $\theta$ and $\phi$, the unit direction of input ray is $(\sin\theta \cos\phi,\sin\theta \sin\phi, \cos\phi )$. For the transmitted ray, $\phi^{\prime} = \phi + \pi$. Applying Snell's law
$$
\sin\theta^{\prime} = \eta \sin\theta
$$
 $\eta$ is the ratio of old index of refraction to the new index of refraction. By Trigonometric Identity Theorem:
$$
\cos\theta ^{\prime} = \sqrt{1-\sin^2\theta^{\prime}} = \sqrt{1-\eta^2 \sin^2\theta} =  \sqrt{1-\eta^2 (1-\cos^2\theta)}
$$
The new $\cos\theta ^{\prime}$ has a chance to be negative, which is a total internal reflection:

1. The light ray must attempt to pass from a medium with a higher refractive index to a medium with a lower refractive index.
2. The angle of incidence of the light ray must be greater than the critical angle.





### Glass

When ray hit a glass surface, both reflection and refraction can occur. As mentioned above, when the $\cos\theta ^{\prime}$ is negative, the total internal reflection occur. Otherwise, using the the Schlick's approximation to determine the probability of the ratio of reflection in pure refraction.
$$
r = \frac{(\eta-1)^2}{(\eta+1)^2}\\
R =r+(1-r)(1-\cos\theta)^5\\
$$
When  $\cos\theta ^{\prime}$ is positive, with a probability of $R$, the ray reflect, and with a probability of $1-R$, the ray refract. Here, $r$ is a intermediate term to calculate $R$.

| Naive Approach                                               | Schlick's Approximation                                      |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| ![image-20230726113141940](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230726113141940.png) | ![image-20230726111049381](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230726111049381.png) |
| ![image-20230726112354899](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230726112354899.png) | ![image-20230726111511208](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230726111511208.png) |

## Depth of Field

In the previous implementation, the model uses a pin-hole camera, which is a sharp focus on everything. To imitate the thin lens, we generate the ray from the camera with a random sample according to the lens property.

1. Rays from the same point on the plane of focus will always be focused to the same point on the image plane, no matter where they pass through the lens.
2. The ray passing though the lens' center won't change direction.

 First, uniformly random generate a point on the lens plane, and get the vector 





| focal distance | Pin-hole                                                     | Thin-Lens                                                    |
| -------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 4.5            | ![image-20230726111049381](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230726111049381.png) | <img src="C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230726144407585.png" alt="image-20230726144407585" style="zoom:150%;" /> |
| 2.2            | ![image-20230726143748988](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230726143748988.png) | ![image-20230726142629082](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230726142629082.png) |













​	









| Spot Light with fog effect                                   | Spot Light without fog effect                                |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| ![image-20230723160707495](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230723160707495.png) | ![image-20230723164316372](C:\Users\surface\AppData\Roaming\Typora\typora-user-images\image-20230723164316372.png) |

















