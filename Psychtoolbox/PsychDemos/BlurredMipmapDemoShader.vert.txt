/* BlurredMipmapDemoShader.vert.txt
 * Shader used for BlurredMipmapDemo.m
 *
 * This shader computes the level of resolution for each output texel
 * during Screen('DrawTexture') of the video image texture. Resolution
 * directly depends on radial distance to the provided simulated center of gaze.
 *
 * Resolution level is used to determine the mip-map miplevel (lod) to use for
 * lookup of the mip-map filtered texel in the images mipmap pyramid.
 *
 * (C) 2012 Mario Kleiner - Licensed under MIT license.
 *
 */

/* Input from Screen('DrawTexture'): */
attribute vec4 auxParameters0;

/* Passed to fragment shader: */
varying vec2  gazePosition;
varying float gazeRadius;
varying vec4  baseColor;

void main(void)
{
    /* Apply standard geometric transformations: */
    gl_Position = ftransform();

    /* Pass standard texture coordinates: */
    gl_TexCoord[0] = gl_MultiTexCoord0;

    /* Pass 'gazePosition' from first two auxParameters: */
    gazePosition.xy = auxParameters0.xy;

    /* Pass 'gazeRadius' from third auxParameters element: */
    gazeRadius = auxParameters0[2];

    /* Base color: */
    baseColor = gl_Color;

    return;
}
