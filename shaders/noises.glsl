#ifndef RAYMARCH_NOISES
#define RAYMARCH_NOISES

#if 1
    // for integer stepped ranges, (ie value-noise/perlin noise functions)
    #define HASHSCALE1 0.1031
    #define HASHSCALE3 vec3(0.1031, 0.1030, 0.0973)
    #define HASHSCALE4 vec4(0.1031, 0.1030, 0.0973, 0.1099)
#else
    // for smaller input rangers (like audio tick or 0-1 UVs use these...)
    #define HASHSCALE1 443.8975
    #define HASHSCALE3 vec3(443.897, 441.423, 437.195)
    #define HASHSCALE4 vec4(443.897, 441.423, 437.195, 444.129)
#endif
float hash11(float p) {
	vec3 p3 = fract(vec3(p) * HASHSCALE1);
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.x+p3.y)*p3.z);
}
float hash12(vec2 p) {
	vec3 p3 = fract(vec3(p.xyx) * HASHSCALE1);
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.x+p3.y)*p3.z);
}
float hash13(vec3 p3) {
	p3 = fract(p3 * HASHSCALE1);
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.x+p3.y)*p3.z);
}
vec2 hash21(float p) {
	vec3 p3 = fract(vec3(p) * HASHSCALE3);
	p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.xx+p3.yz)*p3.zy);
}
vec2 hash22(vec2 p) {
	vec3 p3 = fract(vec3(p.xyx) * HASHSCALE3);
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.xx+p3.yz)*p3.zy);
}
vec2 hash23(vec3 p3) {
	p3 = fract(p3 * HASHSCALE3);
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.xx+p3.yz)*p3.zy);
}
vec3 hash31(float p) {
   vec3 p3 = fract(vec3(p) * HASHSCALE3);
   p3 += dot(p3, p3.yzx+19.19);
   return fract((p3.xxy+p3.yzz)*p3.zyx); 
}
vec3 hash32(vec2 p) {
	vec3 p3 = fract(vec3(p.xyx) * HASHSCALE3);
    p3 += dot(p3, p3.yxz+19.19);
    return fract((p3.xxy+p3.yzz)*p3.zyx);
}
vec3 hash33(vec3 p3) {
	p3 = fract(p3 * HASHSCALE3);
    p3 += dot(p3, p3.yxz+19.19);
    return fract((p3.xxy + p3.yxx)*p3.zyx);
}
vec4 hash41(float p) {
	vec4 p4 = fract(vec4(p) * HASHSCALE4);
    p4 += dot(p4, p4.wzxy+19.19);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);
}
vec4 hash42(vec2 p) {
	vec4 p4 = fract(vec4(p.xyxy) * HASHSCALE4);
    p4 += dot(p4, p4.wzxy+19.19);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);
}
vec4 hash43(vec3 p) {
	vec4 p4 = fract(vec4(p.xyzx)  * HASHSCALE4);
    p4 += dot(p4, p4.wzxy+19.19);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);
}
vec4 hash44(vec4 p4) {
	p4 = fract(p4 * HASHSCALE4);
    p4 += dot(p4, p4.wzxy+19.19);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);
}

//==============================================================================

float noise11(float x) {
    float p = floor(x);
    float f = fract(x);
    f = f*f*(3.0-2.0*f);
    return mix( hash11(p + 0.0), hash11(p + 1.0), f );
}
float noise12(vec2 x) {
    vec2 i = floor(x);
    vec2 f = fract(x);
	float a = hash12(i);
    float b = hash12(i + vec2(1.0, 0.0));
    float c = hash12(i + vec2(0.0, 1.0));
    float d = hash12(i + vec2(1.0, 1.0));
    vec2 u = f*f*(3.0-2.0*f);
	return mix(a, b, u.x) + (c - a) * u.y * (1.0 - u.x) + (d - b) * u.x * u.y;
}
float noise12(sampler2D tex, vec2 x) {
    vec2 p = floor(x);
    vec2 f = fract(x);
	f = f*f*(3.0-2.0*f);
	vec2 uv = p.xy + f.xy;
	return textureLod(tex, (uv + 0.5)/256.0, 0.0).x;
}
float noise13(vec3 x) {
    const vec3 step = vec3(110.0, 241.0, 171.0);
    vec3 i = floor(x);
    vec3 f = fract(x); 
    float n = dot(i, step);
    vec3 u = f*f*(3.0-2.0*f);
    return mix(mix(mix( hash11(n + dot(step, vec3(0.0, 0.0, 0.0))), hash11(n + dot(step, vec3(1.0, 0.0, 0.0))), u.x),
                   mix( hash11(n + dot(step, vec3(0.0, 1.0, 0.0))), hash11(n + dot(step, vec3(1.0, 1.0, 0.0))), u.x), u.y),
               mix(mix( hash11(n + dot(step, vec3(0.0, 0.0, 1.0))), hash11(n + dot(step, vec3(1.0, 0.0, 1.0))), u.x),
                   mix( hash11(n + dot(step, vec3(0.0, 1.0, 1.0))), hash11(n + dot(step, vec3(1.0, 1.0, 1.0))), u.x), u.y), u.z);
}
float noise13(sampler2D tex, vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
	f = f*f*(3.0-2.0*f);
	vec2 uv = (p.xy + vec2(37.0, 17.0)*p.z) + f.xy;
	vec2 rg = textureLod(tex, (uv + 0.5)/256.0, 0.0).yx;
	return mix( rg.x, rg.y, f.z );
}
vec2 noise22(vec2 x) {
    return vec2( noise12(x), noise12(x+17.0) );//OK
    //return vec2( noise12(x), noise12(x-43.0) );//OK
}

//==============================================================================

float Hash11( float p ) {
    return fract( p*17.0*fract( p*0.3183099 ) );
}
float Hash12(vec2 p) {
#if 1
    float h = dot( p, vec2(127.1,311.7) );
    return fract( sin(h) * 43758.5453123 );
#else
    p = 50.0*fract( p*0.3183099 );
    return fract( p.x*p.y*(p.x+p.y) );
#endif
}
float Hash13(vec3 p) {
    // replace this by something better
    p  = 50.0*fract( p*0.3183099 + vec3(0.71,0.113,0.419));
    return fract( p.x*p.y*p.z*(p.x+p.y+p.z) );
}
vec2 Hash21( float p ) {
    return fract( sin(vec2(p, p+1.0)) * vec2(43758.5453123, 22578.1459123) );
}
vec2 Hash22(vec2 p) {
    // replace this by something better
    const vec2 k = vec2( 0.3183099, 0.3678794 );
    p = p*k + k.yx;
    return fract( 16.0 * k*fract( p.x*p.y*(p.x + p.y)) );
}
vec3 Hash33(vec3 p) {
    // replace this by something better
	p = vec3( dot(p,vec3(127.1,311.7, 74.7)),
			  dot(p,vec3(269.5,183.3,246.1)),
			  dot(p,vec3(113.5,271.9,124.6)));
    return fract(sin(p)*43758.5453123);
}

//==============================================================================

float noise(vec2 x)
// value noise 2D
{
#if 0
    vec2 p = floor(x);
    vec2 f = fract(x);
	vec2 u = f*f*(3.0-2.0*f);
    return mix( mix( Hash12( p + vec2(0.0,0.0) ), 
                     Hash12( p + vec2(1.0,0.0) ), u.x),
                mix( Hash12( p + vec2(0.0,1.0) ), 
                     Hash12( p + vec2(1.0,1.0) ), u.x), u.y);
#else
    vec2 p = floor(x);
    vec2 w = fract(x);
    vec2 u = w*w*w*(w*(w*6.0-15.0)+10.0);
    float a = Hash12(p+vec2(0.0,0.0));
    float b = Hash12(p+vec2(1.0,0.0));
    float c = Hash12(p+vec2(0.0,1.0));
    float d = Hash12(p+vec2(1.0,1.0));
    return a + (b-a)*u.x + (c-a)*u.y + (a-b-c+d)*u.x*u.y;
#endif
}
float noise(sampler2D tex, vec2 x)
// value noise 2D
{
    // tex = 256 x 256
    vec2 p = floor(x);
    vec2 f = fract(x);
	f = f*f*(3.0-2.0*f);
#if 0 // low qaulity
	vec2 uv = (p.xy + vec2(37.0,17.0)) + f.xy;
	vec2 rg = textureLod( tex, (uv + 0.5)/256.0, 0.0 ).yx;
#else // high quality
	vec2 uv = (p.xy+vec2(37.0,17.0));
	vec2 rg1 = textureLod( tex, (uv + vec2(0.5,0.5))/256.0, 0.0 ).yx;
	vec2 rg2 = textureLod( tex, (uv + vec2(1.5,0.5))/256.0, 0.0 ).yx;
	vec2 rg3 = textureLod( tex, (uv + vec2(0.5,1.5))/256.0, 0.0 ).yx;
	vec2 rg4 = textureLod( tex, (uv + vec2(1.5,1.5))/256.0, 0.0 ).yx;
	vec2 rg = mix( mix(rg1,rg2,f.x), mix(rg3,rg4,f.x), f.y );
#endif
	return mix( rg.x, rg.y, f.x );
}
float noise(vec3 x)
// value noise 3D
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
    return mix(mix(mix( Hash13( p + vec3(0,0,0) ), 
                        Hash13( p + vec3(1,0,0) ), f.x),
                   mix( Hash13( p + vec3(0,1,0) ), 
                        Hash13( p + vec3(1,1,0) ), f.x), f.y),
               mix(mix( Hash13( p + vec3(0,0,1) ), 
                        Hash13( p + vec3(1,0,1) ), f.x),
                   mix( Hash13( p + vec3(0,1,1) ), 
                        Hash13( p + vec3(1,1,1) ), f.x), f.y), f.z);
}
float noise(sampler2D tex, vec3 x)
// value noise 3D : much faster than the above "float noise(vec3)"
{
	// tex = 256 x 256
    vec3 p = floor(x);
    vec3 f = fract(x);
	f = f*f*(3.0-2.0*f);
#if 1 // low qaulity
	vec2 uv = (p.xy + vec2(37.0,17.0)*p.z) + f.xy;
	vec2 rg = textureLod( tex, (uv + 0.5)/256.0, 0.0 ).yx;
#else // high quality
	vec2 uv = (p.xy+vec2(37.0,17.0)*p.z);
	vec2 rg1 = textureLod( tex, (uv + vec2(0.5,0.5))/256.0, 0.0 ).yx;
	vec2 rg2 = textureLod( tex, (uv + vec2(1.5,0.5))/256.0, 0.0 ).yx;
	vec2 rg3 = textureLod( tex, (uv + vec2(0.5,1.5))/256.0, 0.0 ).yx;
	vec2 rg4 = textureLod( tex, (uv + vec2(1.5,1.5))/256.0, 0.0 ).yx;
	vec2 rg = mix( mix(rg1,rg2,f.x), mix(rg3,rg4,f.x), f.y );
#endif
	return mix( rg.x, rg.y, f.z );
}

float gnoise(vec2 p)
// gradient noise 2D
{
    vec2 i = floor( p );
    vec2 f = fract( p );
	vec2 u = f*f*(3.0-2.0*f);
    return mix( mix( dot( Hash22( i + vec2(0.0,0.0) ), f - vec2(0.0,0.0) ), 
                     dot( Hash22( i + vec2(1.0,0.0) ), f - vec2(1.0,0.0) ), u.x),
                mix( dot( Hash22( i + vec2(0.0,1.0) ), f - vec2(0.0,1.0) ), 
                     dot( Hash22( i + vec2(1.0,1.0) ), f - vec2(1.0,1.0) ), u.x), u.y);
}
float gnoise(vec3 p)
// gradient noise 3D
{
    vec3 i = floor( p );
    vec3 f = fract( p );
	vec3 u = f*f*(3.0-2.0*f);
    return mix( mix( mix( dot( Hash33( i + vec3(0.0,0.0,0.0) ), f - vec3(0.0,0.0,0.0) ), 
                          dot( Hash33( i + vec3(1.0,0.0,0.0) ), f - vec3(1.0,0.0,0.0) ), u.x),
                     mix( dot( Hash33( i + vec3(0.0,1.0,0.0) ), f - vec3(0.0,1.0,0.0) ), 
                          dot( Hash33( i + vec3(1.0,1.0,0.0) ), f - vec3(1.0,1.0,0.0) ), u.x), u.y),
                mix( mix( dot( Hash33( i + vec3(0.0,0.0,1.0) ), f - vec3(0.0,0.0,1.0) ), 
                          dot( Hash33( i + vec3(1.0,0.0,1.0) ), f - vec3(1.0,0.0,1.0) ), u.x),
                     mix( dot( Hash33( i + vec3(0.0,1.0,1.0) ), f - vec3(0.0,1.0,1.0) ), 
                          dot( Hash33( i + vec3(1.0,1.0,1.0) ), f - vec3(1.0,1.0,1.0) ), u.x), u.y), u.z );
}

vec3 noised(vec2 x)
// value noise 2D, derivatives
// return value noise (in x) and its derivatives (in yz)
{
    vec2 p = floor(x);
    vec2 f = fract(x);
#if 1
    // quintic interpolation
    vec2 u = f*f*f*(f*(f*6.0-15.0)+10.0);
    vec2 du = 30.0*f*f*(f*(f-2.0)+1.0);
#else
    // cubic interpolation
    vec2 u = f*f*(3.0-2.0*f);
    vec2 du = 6.0*f*(1.0-f);
#endif
    float a = Hash12( p + vec2(0.0,0.0) );
    float b = Hash12( p + vec2(1.0,0.0) );
    float c = Hash12( p + vec2(0.0,1.0) );
    float d = Hash12( p + vec2(1.0,1.0) );
    float k0 = a;
    float k1 = b - a;
    float k2 = c - a;
    float k4 = a - b - c + d;
    return vec3( k0 + k1*u.x + k2*u.y + k4*u.x*u.y,     // value
                 du * vec2(k1 + k4*u.y, k2 + k4*u.x) ); // derivative
}
vec3 noised(sampler2D tex, vec2 x)
// value noise 2D, derivatives
// return value noise (in x) and its derivatives (in yz)
{
    // tex = 256 X 256
    vec2 f = fract(x);
    vec2 p = floor(x);
#if 0
    // quintic interpolation
    vec2 u = f*f*f*(f*(f*6.0-15.0)+10.0);
    vec2 du = 30.0*f*f*(f*(f-2.0)+1.0);
#else
    // cubic interpolation
    vec2 u = f*f*(3.0-2.0*f);
    vec2 du = 6.0*f*(1.0-f);
#endif
	float a = textureLod( tex, (p+vec2(0.5,0.5))/256.0, 0.0 ).x;
	float b = textureLod( tex, (p+vec2(1.5,0.5))/256.0, 0.0 ).x;
	float c = textureLod( tex, (p+vec2(0.5,1.5))/256.0, 0.0 ).x;
	float d = textureLod( tex, (p+vec2(1.5,1.5))/256.0, 0.0 ).x;
	return vec3(a+(b-a)*u.x+(c-a)*u.y+(a-b-c+d)*u.x*u.y, du*(vec2(b-a,c-a)+(a-b-c+d)*u.yx));
}
vec4 noised(vec3 x)
// value noise 3D (0.0 ~ 1.0), derivatives
// return value noise (in x) and its derivatives (in yzw)
{
    vec3 p = floor(x);
    vec3 w = fract(x);
#if 1
    // quintic interpolation
    vec3 u = w*w*w*(w*(w*6.0-15.0)+10.0);
    vec3 du = 30.0*w*w*(w*(w-2.0)+1.0);
#else
    // cubic interpolation
    vec3 u = w*w*(3.0-2.0*w);
    vec3 du = 6.0*w*(1.0-w);
#endif

#if 0
    float a = Hash13(p+vec3(0.0,0.0,0.0));
    float b = Hash13(p+vec3(1.0,0.0,0.0));
    float c = Hash13(p+vec3(0.0,1.0,0.0));
    float d = Hash13(p+vec3(1.0,1.0,0.0));
    float e = Hash13(p+vec3(0.0,0.0,1.0));
	float f = Hash13(p+vec3(1.0,0.0,1.0));
    float g = Hash13(p+vec3(0.0,1.0,1.0));
    float h = Hash13(p+vec3(1.0,1.0,1.0));
#else
    float n = p.x + 317.0*p.y + 157.0*p.z;
    float a = Hash11(n + 0.0);
    float b = Hash11(n + 1.0);
    float c = Hash11(n + 317.0);
    float d = Hash11(n + 318.0);
    float e = Hash11(n + 157.0);
	float f = Hash11(n + 158.0);
    float g = Hash11(n + 474.0);
    float h = Hash11(n + 475.0);
#endif

    float k0 =   a;
    float k1 =   b - a;
    float k2 =   c - a;
    float k3 =   e - a;
    float k4 =   a - b - c + d;
    float k5 =   a - c - e + g;
    float k6 =   a - b - e + f;
    float k7 = - a + b + c - d + e - f - g + h;
    return vec4( k0 + k1*u.x + k2*u.y + k3*u.z + k4*u.x*u.y + k5*u.y*u.z + k6*u.z*u.x + k7*u.x*u.y*u.z, 
                 du * vec3( k1 + k4*u.y + k6*u.z + k7*u.y*u.z,
                            k2 + k5*u.z + k4*u.x + k7*u.z*u.x,
                            k3 + k6*u.x + k5*u.y + k7*u.x*u.y ) );
}

vec3 gnoised(vec2 p)
// gradient noise 2D, derivatives
// return gradient noise (in x) and its derivatives (in yz)
{
    vec2 i = floor( p );
    vec2 f = fract( p );
#if 1
    // quintic interpolation
    vec2 u = f*f*f*(f*(f*6.0-15.0)+10.0);
    vec2 du = 30.0*f*f*(f*(f-2.0)+1.0);
#else
    // cubic interpolation
    vec2 u = f*f*(3.0-2.0*f);
    vec2 du = 6.0*f*(1.0-f);
#endif
    vec2 ga = Hash22( i + vec2(0.0,0.0) );
    vec2 gb = Hash22( i + vec2(1.0,0.0) );
    vec2 gc = Hash22( i + vec2(0.0,1.0) );
    vec2 gd = Hash22( i + vec2(1.0,1.0) );
    float va = dot( ga, f - vec2(0.0,0.0) );
    float vb = dot( gb, f - vec2(1.0,0.0) );
    float vc = dot( gc, f - vec2(0.0,1.0) );
    float vd = dot( gd, f - vec2(1.0,1.0) );
    return vec3( va + u.x*(vb-va) + u.y*(vc-va) + u.x*u.y*(va-vb-vc+vd),   // value
                 ga + u.x*(gb-ga) + u.y*(gc-ga) + u.x*u.y*(ga-gb-gc+gd) +  // derivatives
                 du * (u.yx*(va-vb-vc+vd) + vec2(vb,vc) - va));
}
vec4 gnoised(vec3 x)
// gradient noise 3D, derivatives
// return value noise (in x) and its derivatives (in yzw)
{
    // grid
    vec3 p = floor(x);
    vec3 w = fract(x);
    #if 1
    // quintic interpolant
    vec3 u = w*w*w*(w*(w*6.0-15.0)+10.0);
    vec3 du = 30.0*w*w*(w*(w-2.0)+1.0);
    #else
    // cubic interpolant
    vec3 u = w*w*(3.0-2.0*w);
    vec3 du = 6.0*w*(1.0-w);
    #endif
    // gradients
    vec3 ga = Hash33( p+vec3(0.0,0.0,0.0) );
    vec3 gb = Hash33( p+vec3(1.0,0.0,0.0) );
    vec3 gc = Hash33( p+vec3(0.0,1.0,0.0) );
    vec3 gd = Hash33( p+vec3(1.0,1.0,0.0) );
    vec3 ge = Hash33( p+vec3(0.0,0.0,1.0) );
	vec3 gf = Hash33( p+vec3(1.0,0.0,1.0) );
    vec3 gg = Hash33( p+vec3(0.0,1.0,1.0) );
    vec3 gh = Hash33( p+vec3(1.0,1.0,1.0) );
    // projections
    float va = dot( ga, w-vec3(0.0,0.0,0.0) );
    float vb = dot( gb, w-vec3(1.0,0.0,0.0) );
    float vc = dot( gc, w-vec3(0.0,1.0,0.0) );
    float vd = dot( gd, w-vec3(1.0,1.0,0.0) );
    float ve = dot( ge, w-vec3(0.0,0.0,1.0) );
    float vf = dot( gf, w-vec3(1.0,0.0,1.0) );
    float vg = dot( gg, w-vec3(0.0,1.0,1.0) );
    float vh = dot( gh, w-vec3(1.0,1.0,1.0) );
    // interpolations
    return vec4( va + u.x*(vb-va) + u.y*(vc-va) + u.z*(ve-va) + u.x*u.y*(va-vb-vc+vd) + u.y*u.z*(va-vc-ve+vg) + u.z*u.x*(va-vb-ve+vf) + (-va+vb+vc-vd+ve-vf-vg+vh)*u.x*u.y*u.z,    // value
                 ga + u.x*(gb-ga) + u.y*(gc-ga) + u.z*(ge-ga) + u.x*u.y*(ga-gb-gc+gd) + u.y*u.z*(ga-gc-ge+gg) + u.z*u.x*(ga-gb-ge+gf) + (-ga+gb+gc-gd+ge-gf-gg+gh)*u.x*u.y*u.z +   // derivatives
                 du * (vec3(vb,vc,ve) - va + u.yzx*vec3(va-vb-vc+vd,va-vc-ve+vg,va-vb-ve+vf) + u.zxy*vec3(va-vb-ve+vf,va-vb-vc+vd,va-vc-ve+vg) + u.yzx*u.zxy*(-va+vb+vc-vd+ve-vf-vg+vh) ));
}

//==============================================================================

#define NUM_NOISE_OCTAVES   5

float fbm11(float x) {
	float v = 0.0;
	float a = 0.5;
	float shift = float(100.0);
	for (int i = 0; i < NUM_NOISE_OCTAVES; i++) {
		v += a * noise11(x);
		x = x * 2.0 + shift;
		a *= 0.5;
	}
	return v;
}
float fbm12(vec2 x) {
	float v = 0.0;
	float a = 0.5;
	vec2 shift = vec2(100.0);
    mat2 rot = mat2(cos(0.5), sin(0.5), -sin(0.5), cos(0.5));
	for (int i = 0; i < NUM_NOISE_OCTAVES; i++) {
		v += a * noise12(x);
		x = rot * x * 2.0 + shift;
		a *= 0.5;
	}
	return v;
}
float fbm13(vec3 x) {
	float v = 0.0;
	float a = 0.5;
	vec3 shift = vec3(100.0);
	for (int i = 0; i < NUM_NOISE_OCTAVES; i++) {
		v += a * noise13(x);
		x = x * 2.0 + shift;
		a *= 0.5;
	}
	return v;
}

const mat2 rotM2  = mat2(0.80,0.60,-0.60,0.80);
const mat2 rotM2i = mat2(0.80,-0.60,0.60,0.80);//36.87(deg)
const mat3 rotM3  = mat3(0.00,0.80,0.60,-0.80,0.36,-0.48,-0.60,-0.48,0.64);
const mat3 rotM3i = mat3(0.00,-0.80,-0.60,0.80,0.36,-0.48,0.60,-0.48,0.64);

float fbm(sampler2D tex, vec2 p) {
    // tex = 256 x 256 (grayNoise256.png)
    float f = 0.0;
    f += 0.5000 * texture( tex, p/256.0 ).x; p = rotM2*p*2.02;
    f += 0.2500 * texture( tex, p/256.0 ).x; p = rotM2*p*2.03;
    f += 0.1250 * texture( tex, p/256.0 ).x; p = rotM2*p*2.01;
    f += 0.0625 * texture( tex, p/256.0 ).x;
    return f/0.9375;
}

float fbm_4(vec2 x) {
    float f = 1.92;
    float a = 0.0;
    float b = 0.5;
    for(int i=0; i<4; i++)
    {
        float n = noise(x);
        a += b * n;
        b *= 0.5;//0.55;
        x = f * rotM2 * x;
    }
	return a;
}
float fbm_4(sampler2D tex, vec2 x) {
    float f = 1.92;//2.0;
    float a = 0.0;
    float b = 0.5;
    for(int i=0; i<4; i++)
    {
        //float n = noise12(tex, x); // low quality
        float n = noise(tex, x);     // high quality
        a += b * n;
        b *= 0.5;//0.55;
        x = f * rotM2 * x;
    }
	return a;
}
float fbm_9(vec2 x) {
    float f = 1.92;//2.0;
    float a = 0.0;
    float b = 0.5;
    for(int i=0; i<9; i++)
    {
        float n = noise(x);
        a += b * n;
        b *= 0.5;//0.55;
        x = f * rotM2 * x;
    }
	return a;
}
float fbm_9(sampler2D tex, vec2 x) {
    float f = 1.92;//2.0;
    float a = 0.0;
    float b = 0.5;
    for(int i=0; i<9; i++)
    {
        //float n = noise12(tex, x); //low quality
        float n = noise(tex, x);     //high quality
        a += b * n;
        b *= 0.5;//0.55;
        x = f * rotM2 * x;
    }
	return a;
}

float fbm_4(vec3 x) {
    float f = 1.92;//2.0;
    float a = 0.0;
    float b = 0.5;
    for(int i=0; i<4; i++ )
    {
        float n = noise(x);
        a += b * n;
        b *= 0.5;
        x = f * rotM3 * x;
    }
	return a;
}
float fbm_4(sampler2D tex, vec3 x) {
    float f = 1.92;//2.0;
    float a = 0.0;
    float b = 0.5;
    for(int i=0; i<4; i++ )
    {
        float n = noise(tex, x);
        a += b * n;
        b *= 0.5;
        x = f * rotM3 * x;
    }
	return a;
}

vec3 fbmd_5(vec2 x) {
    float f = 1.92;
    float a = 0.0;
    float b = 0.5;
    vec2  d = vec2(0.0);
    mat2  m = mat2(1.0,0.0, 0.0,1.0);
    for(int i=0; i<5; i++)
    {
        vec3 n = noised(x);
        a += b * n.x;          // accumulate values		
        d += b * m * n.yz;       // accumulate derivatives
        b *= 0.5;
        x = f * rotM2 * x;
        m = f * rotM2i * m;
    }
	return vec3( a, d );
}
vec3 fbmd_9(vec2 x) {
    float f = 1.92;
    float a = 0.0;
    float b = 0.5;
    vec2  d = vec2(0.0);
    mat2  m = mat2(1.0,0.0, 0.0,1.0);
    for(int i=0; i<9; i++)
    {
        vec3 n = noised(x);
        a += b * n.x;          // accumulate values		
        d += b * m * n.yz;       // accumulate derivatives
        b *= 0.5;
        x = f * rotM2 * x;
        m = f * rotM2i * m;
    }
	return vec3( a, d );
}

vec4 fbmd_5(vec3 x) {
    // return value (in x) and its derivatives (in yzw)
    float f = 1.92;
    float a = 0.0;
    float b = 0.5;
    vec3  d = vec3(0.0);
    mat3  m = mat3(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0);
    for(int i=0; i<5; i++)
    {
        vec4 n = noised(x);
        a += b * n.x;       // accumulate values		
        d += b * m * n.yzw; // accumulate derivatives
        b *= 0.5;
        x = f * rotM3 * x;
        m = f * rotM3i * m;
    }
	return vec4( a, d );
}
vec4 fbmd_7(vec3 x) {
    // return value (in x) and its derivatives (in yzw)
    float f = 1.92;
    float a = 0.0;
    float b = 0.5;
    vec3  d = vec3(0.0);
    mat3  m = mat3(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0);
    for(int i=0; i<7; i++)
    {
        vec4 n = noised(x);
        a += b * n.x;       // accumulate values		
        d += b * m * n.yzw; // accumulate derivatives
        b *= 0.5;
        x = f * rotM3 * x;
        m = f * rotM3i * m;
    }
	return vec4( a, d );
}

//==============================================================================

float triSmooth11(float x){ return 0.25+0.25*cos(PI2*x); }
float triSmooth12(vec2 p) { return triSmooth11(p.x+triSmooth11(p.y)) + triSmooth11(p.y+triSmooth11(p.x)); }
vec3  triSmooth33(vec3 p) { return vec3(triSmooth11(p.x+triSmooth11(p.y)), triSmooth11(p.y+triSmooth11(p.z)), triSmooth11(p.z+triSmooth11(p.x))); }
float triNoiseSmooth13(vec3 p) {
    float z = 1.4;
	float rz = 0.0;
    vec3 bp = p;
	for(int i = 0; i <= 3; i++)
	{
        vec3 dg = triSmooth33( bp );
        p += dg;
        bp *= 1.8;//2.0
		z *= 1.5;
		p *= 1.2;
        rz += (triSmooth11(p.z + triSmooth11(p.x + triSmooth11(p.y))))/z;
        bp += 0.14;
	}
	return rz;
}
float tri11(float x) { return abs(fract(x)-0.5); }
float tri12(vec2 p)  { return tri11(p.x+tri11(p.y)) + tri11(p.y+tri11(p.x)); }
vec3  tri33(vec3 p)  { return vec3(tri11(p.x+tri11(p.y)), tri11(p.y+tri11(p.z)), tri11(p.z+tri11(p.x))); }
float triNoise13(vec3 p) {
    float z = 1.4;
	float rz = 0.0;
    vec3 bp = p;
	for(int i = 0; i <= 3; i++)
	{
        vec3 dg = tri33( bp );
        p += dg;
        bp *= 1.8;//2.0
		z *= 1.5;
		p *= 1.2;
        rz += (tri11(p.z + tri11(p.x + tri11(p.y))))/z;
        bp += 0.14;
	}
	return rz;
}

//==============================================================================

float sinusoidBumps(in vec3 p, in float time) {
    // NOTE*: sinusoidBumps() + combined with 3D rotations ==> better than cheap noise
    return sin(p.x*16.+time*0.57)*cos(p.y*16.+time*2.17)*sin(p.z*16.-time*1.31) + 0.5*sin(p.x*32.+time*0.07)*cos(p.y*32.+time*2.11)*sin(p.z*32.-time*1.23);
}

#endif // RAYMARCH_NOISES