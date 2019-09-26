#ifndef RAYMARCH_DISTANCES
#define RAYMARCH_DISTANCES

//==============================================================================
// distance field functions (primitives)
//==============================================================================

float fPlane(vec3 p) { return p.y; }
float fSphere(vec3 p, float s) { return length(p)-s; }
float fBox(vec3 p, vec3 b) {
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}
float fEllipsoid(vec3 p, vec3 r) { return (length( p/r ) - 1.0) * min(min(r.x,r.y),r.z); }
float uRoundBox(vec3 p, vec3 b, float r) { return length(max(abs(p)-b,0.0))-r; }
float fRoundBox(vec3 p, vec3 b, float r) { return length(max(abs(p)-b+vec3(r),0.0))-r; }
float fTorus(vec3 p, vec2 t) { return length( vec2(length(p.xz)-t.x,p.y) )-t.y; }
float fHexPrism(vec3 p, vec2 h) {
    vec3 q = abs(p);
#if 0
    return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
#else
    float d1 = q.z-h.y;
    float d2 = max((q.x*0.866025+q.y*0.5),q.y)-h.x;
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
#endif
}
float fCapsule(vec3 p, vec3 a, vec3 b, float r) {
	vec3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return length( pa - ba*h ) - r;
}
float fEquilateralTriangle(vec2 p) {
    const float k = 1.73205;//sqrt(3.0);
    p.x = abs(p.x) - 1.0;
    p.y = p.y + 1.0/k;
    if( p.x + k*p.y > 0.0 ) p = vec2( p.x - k*p.y, -k*p.x - p.y )/2.0;
    p.x += 2.0 - 2.0*clamp( (p.x+2.0)/2.0, 0.0, 1.0 );
    return -length(p)*sign(p.y);
}
float fTriPrism(vec3 p, vec2 h) {
    vec3 q = abs(p);
    float d1 = q.z-h.y;
#if 1
    // distance bound
    float d2 = max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5;
#else
    // correct distance
    h.x *= 0.866025;
    float d2 = fEquilateralTriangle(p.xy/h.x)*h.x;
#endif
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}
float fCylinder(vec3 p, vec2 h) {
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}
float fCone(vec3 p, vec3 c) {
    vec2 q = vec2( length(p.xz), p.y );
    float d1 = -q.y-c.z;
    float d2 = max( dot(q,c.xy), q.y);
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}
float fConeSection(vec3 p, float h, float r1, float r2) {
    float d1 = -p.y - h;
    float q = p.y - h;
    float si = 0.5*(r1-r2)/h;
    float d2 = max( sqrt( dot(p.xz,p.xz)*(1.0-si*si)) + q*si - r2, q );
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}
float fPryamid4(vec3 p, vec3 h) {
    // h = (cos a, sin a, height)
    // Tetrahedron = Octahedron - Cube
    float box = fBox( p - vec3(0,-2.0*h.z,0), vec3(2.0*h.z) );
    float d = 0.0;
    d = max( d, abs( dot(p, vec3( -h.x, h.y, 0 )) ));
    d = max( d, abs( dot(p, vec3(  h.x, h.y, 0 )) ));
    d = max( d, abs( dot(p, vec3(  0, h.y, h.x )) ));
    d = max( d, abs( dot(p, vec3(  0, h.y,-h.x )) ));
    float octa = d - h.z;
    return max(-box,octa); // Subtraction
}
float length2(vec2 p) { return sqrt( p.x*p.x + p.y*p.y ); }
float length6(vec2 p) {
	p = p*p*p; p = p*p;
	return pow( p.x + p.y, 1.0/6.0 );
}
float length8(vec2 p) {
	p = p*p; p = p*p; p = p*p;
	return pow( p.x + p.y, 1.0/8.0 );
}
float fTorus82(vec3 p, vec2 t) {
    vec2 q = vec2(length2(p.xz)-t.x,p.y);
    return length8(q)-t.y;
}
float fTorus88(vec3 p, vec2 t) {
    vec2 q = vec2(length8(p.xz)-t.x,p.y);
    return length8(q)-t.y;
}
float fCylinder6(vec3 p, vec2 h) { return max( length6(p.xz)-h.x, abs(p.y)-h.y ); }
// float fNoise(vec2 p) {
// 	vec2 s = sin(p * 0.6345) + sin(p * 1.62423);
// 	return dot(s, vec2(0.125)) + 0.5;
// }
float fSinusoidalPlasma(vec3 p) { return sin(p.x)*cos(p.y)*sin(p.z) + 0.5*sin(p.x*2.0)*cos(p.y*2.0)*sin(p.z*2.0); }
float fsmin(float a, float b, float k) {
    // polynomial smooth min (k = 0.1)
	float h = clamp(0.5 + 0.5*(b-a)/k, 0.0, 1.0);
	return mix(b, a, h) - k*h*(1.0-h);
}
#if 1
float fsmax(float a, float b, float k) {
    // polynomial smooth max (k = 0.1)
	float h = clamp(0.5 + 0.5*(b-a)/k, 0.0, 1.0);
	return mix(a, b, h) + k*h*(1.0-h);
}
#else
float fsmax(float a, float b, float k) {
    // polynomial smooth max (k = 0.1)
	float h = clamp(0.5 + 0.5*(a-b)/k, 0.0, 1.0);
	return mix(b, a, h) + k*h*(1.0-h);
    // NOTE: this might be the same result as the above fsmax()
}
#endif

//==============================================================================
// (object-level) operations (vector-valued object: d.xy=(dist, matID))
//==============================================================================

//vec2 dUnion(vec2 d1, vec2 d2) { return d1.x < d2.x ? d1 : d2; }
vec2 dUnion(vec2 d1, vec2 d2) { return mix(d1, d2, step(d2.x, d1.x)); }
vec2 dIntersect(vec2 d1, vec2 d2) { return mix(d2, d1, step(d2.x, d1.x)); }
vec2 dSubtract(vec2 d1, vec2 d2) { return dIntersect(d1, vec2(-d2.x, d2.y)); }
// ds = displacement to apply to d1
// for example:
// float ds(vec3 p) { return sin(2*p.x)*sin(2*p.y)*sin(2*p.z); }
// vec2 new_sphere = dDisplace( vec2(fSphere(p,r), m), ds(p) );
vec2 dDisplace(vec2 d1, float ds) { return vec2(d1.x + ds, d1.y); }
vec2 dSmoothUnion(vec2 d1, vec2 d2, float k) {
    // smooth union (often k = 0.1)
	float f = fsmin( d1.x, d2.x, k );
	float m = mix( d1.y, d2.y, step(d2.x, d1.x) );
	return vec2(f, m);
}
vec2 dSmoothIntersect(vec2 d1, vec2 d2, float k) {
    // smooth intersection (often k = 0.1)
	float f = fsmax( d1.x, d2.x, k );
	float m = mix( d1.y, d2.y, step(d1.x, d2.x) );
	return vec2(f, m);
}
vec2 dSmoothSubtract(vec2 d1, vec2 d2, float k) { return dSmoothIntersect(d1, vec2(-d2.x, d2.y), k); }

//==============================================================================
// (space-level) operations (vector-valued space: p=(x,y,z))
//==============================================================================

vec3 sRepeat(vec3 p, vec3 spacing) {
    // The same domain space is repeated along (x,y,z)-axis, whose size is spacing.xyz, respectively.
    // if spacing.y=0, then infinite column is repeated along x-axis & z-axis
	vec3 q = p - 0.5*spacing;
	return mod(q, spacing) - 0.5*spacing;
}
vec3 sTwist(vec3 p, float angle) {
    // twist by angle per unit-length along y-axis with the right-hand rule
	float c = cos( angle*p.y );
	float s = sin( angle*p.y );
	mat2 m = mat2( c, -s, s, c );
	vec2 twist = m*p.zx;
	return vec3( twist.y, p.y, twist.x );
}
vec3 sTwistY(vec3 p, float angle) {
    // twist along y-axis
	return sTwist( p, angle );
}
vec3 sTwistZ(vec3 p, float angle) {
    // twist along z-axis
	vec3 q = p.yzx;
	q = sTwist( q, angle );
	return q.zxy;
}
vec3 sTwistX(vec3 p, float angle) {
    // twist along x-axis
	vec3 q = p.zxy;
	q = sTwist( q, angle );
	return q.yzx;
}
vec3 sTranslate(vec3 p, vec3 t) { return p - t; }
vec3 sRotateX(vec3 p, float angle) {
	float s = sin(angle);
	float c = cos(angle);
	return vec3( p.x, c*p.y+s*p.z, -s*p.y+c*p.z );
}
vec3 sRotateY(vec3 p, float angle) {
	//float s = sin(angle);
    float s = -sin(angle);
	float c = cos(angle);
	return vec3( c*p.x+s*p.z, p.y, -s*p.x+c*p.z );
}
vec3 sRotateZ(vec3 p, float angle) {
	float s = sin(angle);
	float c = cos(angle);
	return vec3( c*p.x+s*p.y, -s*p.x+c*p.y, p.z );
}
// Scale
// e.g., float fScaled = fSphere( p/s, r ) * s
// only possible to scale uniformly
vec3 sCheapBendY(vec3 p, float angle) {
    // bending moment: yaxis, convex up = zaxis, fixed: (0,0,0)
	float c = cos( angle*p.x );
	float s = sin( angle*p.x );
	mat2 m = mat2( c, -s, s, c );
	vec2 b = m*p.zx;
	return vec3( b.y, p.y, b.x );
}
vec3 sCheapBendZ(vec3 p, float angle) {
    // bending moment: zaxis, convex up = xaxis, fixed: (0,0,0)
	float c = cos( angle*p.y );
	float s = sin( angle*p.y );
	mat2 m = mat2( c, -s, s, c );
	vec2 b = m*p.xy;
	return vec3( b.x, b.y, p.z );
}
vec3 sCheapBendX(vec3 p, float angle) {
    // bending moment: xaxis, convex up: yaxis, fixed: (0,0,0)
	float c = cos( angle*p.z );
	float s = sin( angle*p.z );
	mat2 m = mat2( c, -s, s, c );
	vec2 b = m*p.yz;
	return vec3( p.x, b.x, b.y );
}

//==============================================================================
// Ray marching
//==============================================================================

#define EPS             0.001
#define MAX_RAY_ITER    512     // 256 128
#define MIN_RAY_DIST    0.02    // 0.02
#define MAX_RAY_DIST    20.0    // 20.0
#define MAX_SHADOW_ITER 64      // 16
#define MIN_SHADOW_DIST 0.005   // 0.01
#define MAX_SHADOW_DIST 10.0    // 2.5

vec2 sceneMap( in vec3 p );

//==============================================================================

float sceneShadow( in vec3 ro, in vec3 rd, in float tmin, in float tmax, in float k )
// tmin = 0.005
// tmax = MAX_RAY_DIST*0.5
// k = fade-off factor (eg: 10, 32 or 64) (smaller: soft penumbra, larger: hard-edged penumbra)
{
	float shade = 1.0;
    float t = tmin;
    for(int i = 0; i < MAX_SHADOW_ITER; i++)
    {
        float h = sceneMap(ro + rd * t).x;
        shade = min( shade, k * h / t );
        t += clamp( h, 0.5, 1.0 );
        if( h < EPS || t > tmax ) break;
    }
    return min( max(shade, 0.0) + 0.1, 1.0 ); // 0.3 or 0.1 (preference)
}
float sceneShadow( in vec3 ro, in vec3 rd )
{
    return sceneShadow( ro, rd, MIN_SHADOW_DIST, MAX_RAY_DIST*0.5, 64.0 ); // 10, 32 or 64
}

float sceneShadow( in vec3 ro, in vec3 rd, float tmin, float k )
// k = fade-off factor (eg: 10, 32 or 64) (smaller: soft penumbra, larger: hard-edged penumbra)
{
    float shade = 1.0;
    float t = tmin;
    for(int i=0; i<40; i++)
    {
        float h = sceneMap(ro + rd*t).x;
        shade = min( shade, smoothstep(0.0, 1.0, k*h/t) );
		t += clamp( h, 0.05, 0.5 );
		if( h < 0.0001 ) break;
    }
    return clamp(shade, 0.0, 1.0);
}

//==============================================================================

float sceneAO( in vec3 p, in vec3 n )
{
    float ao = 0.0;
    float s = 1.0;
    for(int i = 0; i < 6; i++)
    {
        float off = 0.001 + 0.2 * float(i)/5.0;
        float t = sceneMap( n * off + p ).x;
        ao += ( off - t ) * s;
        s *= 0.4;
    }
    return smoothstep( 0.0, 1.0, clamp(1.0-12.0*ao, 0.0, 1.0) );
}

//==============================================================================

vec3 sceneNormal( in vec3 p )
{
    vec3 eps = vec3(EPS, 0.0, 0.0);
    vec3 n = vec3(
            sceneMap(p + eps.xyy).x - sceneMap(p - eps.xyy).x,
            sceneMap(p + eps.yxy).x - sceneMap(p - eps.yxy).x,
            sceneMap(p + eps.yyx).x - sceneMap(p - eps.yyx).x);
    return normalize(n);
}

vec3 sceneNormal( in vec3 p, in float t )
{
    vec2 eps = vec2( 0.005*t, 0.0 );
	return normalize( vec3(
            sceneMap(p + eps.xyy).x - sceneMap(p - eps.xyy).x,
            sceneMap(p + eps.yxy).x - sceneMap(p - eps.yxy).x,
            sceneMap(p + eps.yyx).x - sceneMap(p - eps.yyx).x ) );
}

//==============================================================================

float rayTracing( in vec3 ro, in vec3 rd, in float tmin, in float tmax, in float mfac )
{
    float t = tmin;
	for(int i = 0; i < MAX_RAY_ITER; i++)
	{
		float d = sceneMap( ro + rd*t ).x;
		if( d < EPS || t > tmax ) break;
        t += mfac * d; // mfac = 1.0, 0.75, 0.5, 0.3...
	}
	return t;
}

vec2 rayMarching( in vec3 ro, in vec3 rd, float tmin, float tmax )
// return vec2(travelDist, materialID)
// if travelDist > tmax, no intersection found
{
    float m = -1.0;
    float t = tmin;
    for(int i = 0; i < MAX_RAY_ITER; i++)
    {
        vec2 d = sceneMap( ro + rd*t );
        if( d.x < EPS || t > tmax ) break;
        t += d.x;
        m = d.y;
    }
    if( t > tmax ) m = -1.0;
    return vec2(t, m);
}

vec2 rayMarching( in vec3 ro, in vec3 rd )
{
    return rayMarching( ro, rd, MIN_RAY_DIST, MAX_RAY_DIST );
}

//==============================================================================

float bisectTracing( in vec3 ro, in vec3 rd, in float tmin, in float tmax )
{
    float t = tmin, told = tmin, mid, dn;
    float d = sceneMap(ro + rd*t).x;
    float sgn = sign(d);
    for(int i=0; i<80; i++)
    {
        if( sign(d) != sgn || d < EPS || t > tmax ) break;
        told = t;
        t += step(-1.0, -d)*(log(abs(d) + 1.1)*0.7 - d*0.7) + d*0.7;
        d = sceneMap(ro + rd*t).x;
    }
    if( sign(d) != sgn )
    {
        dn = sign( sceneMap(ro + rd*told).x );
        vec2 iv = vec2(told, t);
        for(int ii=0; ii<8; ii++)
        {
            mid = dot(iv, vec2(0.5));
            float d = sceneMap(ro + rd*mid).x;
            if( abs(d) < EPS ) break;
            iv = mix( vec2(iv.x, mid), vec2(mid, iv.y), step(0.0, d*dn) );
        }
        t = mid;
    }
    return t;
}

vec2 bisectMarching( in vec3 ro, in vec3 rd, in float tmin, in float tmax )
// return vec2(travelDist, materialID)
// if travelDist > tmax, no intersection found
{
    float m = -1.0;
    float t = tmin, told = tmin, mid, dn;
    vec2 d = sceneMap(ro + rd*t); m = d.y;
    float sgn = sign(d.x);
    for(int i=0; i<80; i++)
    {
        if( sign(d.x) != sgn || d.x < EPS || t > tmax ) break;
        told = t;
        t += step(-1.0, -d.x)*(log(abs(d.x) + 1.1)*0.7 - d.x*0.7) + d.x*0.7;
        d = sceneMap(ro + rd*t); m = d.y;
    }
    if( sign(d.x) != sgn )
    {
        dn = sign( sceneMap(ro + rd*told).x );
        vec2 iv = vec2(told, t);
        for(int ii=0; ii<8; ii++)
        {
            mid = dot(iv, vec2(0.5));
            vec2 d = sceneMap(ro + rd*mid); m = d.y;
            if( abs(d.x) < EPS ) break;
            iv = mix( vec2(iv.x, mid), vec2(mid, iv.y), step(0.0, d.x*dn) );
        }
        t = mid;
    }
    return vec2(t, m);
}

//==============================================================================

mat3 cameraMatrix( vec3 ro, vec3 ta )
{
    vec3 cw = normalize( ta - ro );
    vec3 cu = normalize( cross(cw, vec3(0.0, 1.0, 0.0)) );
    vec3 cv = normalize( cross(cu, cw) );
    return mat3(cu, cv, cw);
}

#endif // RAYMARCH_DISTANCES