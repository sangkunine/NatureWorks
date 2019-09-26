#ifndef RAYMARCH_RAY
#define RAYMARCH_RAY

// screen position [-1, 1]
vec2 ndc = ( gl_FragCoord.xy * 2.0 - resolution ) / resolution;

// ray direction in normalized device coordinate
vec4 ndcRay = vec4( ndc.xy, 1.0, 1.0 );

// ray direction in world coordinate
ndcRay = cameraProjectionMatrixInverse * ndcRay;
ndcRay = cameraWorldMatrix * ndcRay;
ndcRay /= ndcRay.w;
vec3 rd = normalize( ndcRay.xyz );

// ray origin (= camera position in world coordinate)
vec3 ro = cameraPosition;

#ifdef USE_SHADERTOY_CAMERA
    vec3 ta = (cameraWorldMatrix * cameraProjectionMatrixInverse * vec4(0.0,0.0,0.0,1.0)).xyz;
    vec3 cw = normalize( ta - ro );
    vec3 cv = vec3(0.0, 1.0, 0.0);
    vec3 cu = cross(cw, cv);
    cv = cross(cu, cw);
    mat3 ca = mat3(cu, cv, cw);
#endif

#endif // RAYMARCH_RAY