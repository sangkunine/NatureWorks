#ifndef RAYMARCH_INCLUDE
#define RAYMARCH_INCLUDE

#define gl_FragCoord     (gl_FragCoord / devicePixelRatio)

#define PI 3.14159265359
#define PI2 6.28318530718
#define PI_HALF 1.5707963267949
#define RECIPROCAL_PI 0.31830988618
#define RECIPROCAL_PI2 0.15915494
#define LOG2 1.442695
#define EPSILON 1e-6

#define saturate(a) clamp( a, 0.0, 1.0 )
#define whiteCompliment(a) ( 1.0 - saturate( a ) )

uniform float devicePixelRatio;
uniform float time;
uniform vec2 resolution;
uniform float cameraNear;
uniform float cameraFar;
// uniform mat4 viewMatrix;

// uniform vec3 cameraPosition;
uniform mat4 cameraWorldMatrix;
uniform mat4 cameraProjectionMatrixInverse;
uniform sampler2D textureMaps[4];

float getFragDepth( vec3 fragPosW )
// gl_FragDepthEXT = getFragDepth( fragPos );
// where gl_FragDepthEXT = [0.0, 1.0]
// if gl_FragDepthEXT = 1.0  ==> fragment will be killed
// if gl_FragDepthEXT = 0.99 ==> fragment will be preserved as background
{
    vec4 fragPosV = viewMatrix * vec4(fragPosW, 1.0);

    // // viewZ to perspectiveDepth
    // float fragDepth = cameraFar * (cameraNear + fragPosV.z) / ((cameraFar - cameraNear) * fragPosV.z);

    // viewZ To orthographicDepth
    float fragDepth = ( fragPosV.z + cameraNear ) / ( cameraNear - cameraFar );

    return fragDepth; // orthographicDepth(O), perspectiveDepth(X)
}

#endif // RAYMARCH_INCLUDE