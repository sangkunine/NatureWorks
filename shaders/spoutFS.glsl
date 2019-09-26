@import ./include;

//#define ENABLE_TEST_RENDER

#ifdef LOW_QUALITY
	#define kRaymarchMaxIter 16
#else
	#define kRaymarchMaxIter 32
	#define ENABLE_AMBIENT_OCCLUSION
	#define DOUBLE_SIDED_TRANSPARENCY
#endif

#ifdef ENABLE_HARD_SHADOWS
	#undef ENABLE_SOFT_SHADOWS
#else
	#define ENABLE_SOFT_SHADOWS
#endif

#define ENABLE_SPECULAR
#define ENABLE_REFLECTIONS
#define ENABLE_TRANSPARENCY
#define ENABLE_FOG
#define ENABLE_TEXTURE			//psk: now not working

#define ENABLE_DIRECTIONAL_LIGHT
#define ENABLE_DIRECTIONAL_LIGHT_FLARE
//#define ENABLE_POINT_LIGHT
//#define ENABLE_POINT_LIGHT_FLARE

const float kFarClip = 100.0;
const float kRaymarchEpsilon = 0.001;	//0.01
const float kFogDensity = 0.005;		//0.05
const float kPI = 3.141592654;
const float kTwoPI = kPI * 2.0;
const float kTranspNo = -1.0;
const float kTranspYes = 1.0;
const float kTranspInverse = 0.0;

const float kMaterialGround = 1.0;		// sky=0.0
const float kMaterialGold = 2.0;
const float kMaterialSilver = 3.0;
const float kMaterialWall = 4.0;
const float kMaterialPipe = 5.0;
const float kMaterialWater = 6.0;
const float kMaterialTexture0 = 7.0;	// textureMaps[0]
const float kMaterialTexture1 = 8.0;	// textureMaps[1]

struct CRay
{
	vec3 vOrigin;
	vec3 vDir;
	float fStartDist;
	float fLength;	// maximum distance of ray
};
struct CHitInfo
{
	vec3 vPos;		// vPos = ro + rd * fDist
	float fDist;
	vec3 vObjectID;	// x=ID, yz=??(tex_uv)
}; 
struct CSurface
{
	vec3 vNormal;
	vec3 cReflection;	// reflected color
	vec3 cTransmission;	// transmitted color (refraction)
};
struct CMaterial
{
	vec3 cAlbedo;		// diffuse reflectivity (= cTexture * cTexture)
	float fR0;			// reflection coeff for light incoming parallel to the normal (= [(n1-n2)/(n1+n2)]^2 )
	float fSmoothness;	// related to alpha for (n*h)^alpha
	vec2 vParam;
	float fTransparency;
	float fRefractIndex;	// defined as n1/n2 (where n1=come, n2=go)
	// cf: n(Air) = 1.0003, n(Water)=1.333, n(common Glass)=1.52, n(Diamond)=2.42
	// n(Aluminium/Al)=1.0972, n(Copper/Cu)=0.46090, n(Iron/Fe)=2.9304
	// n(Gold/Au)=0.27049, n(Silver/Ag)=0.15016, n(Titanium/Ti)=2.6112, n(Lead/Pb)=2.6
};
struct CShading
{
	vec3 cDiffuse;
	vec3 cSpecular;
};
struct CPointLight
{
	vec3 vPos;
	vec3 cCol;
};
struct CDirectionalLight
{
	vec3 vDir;
	vec3 cCol;
};

//----------------------------------------------------------------------------
// primitives
//----------------------------------------------------------------------------
float fSphere( vec3 p, float r )
// r = radius
{
	return length(p) - r;
}
float fEllipsoid( in vec3 p, in vec3 r )
{
	return (length( p/r ) - 1.0) * min(min(r.x, r.y), r.z);
}
float fBox( vec3 p, vec3 b )
// b.x =half of width (along x-axis)
// b.y = halft of height (along y-axis)
// b.z = half of depth (along z-axis)
{
	vec3 d = abs(p) - b;
	return min( max(d.x, max(d.y,d.z)), 0.0) + length(max(d,0.0));
}
float fTorus( vec3 p, vec2 t )		// t=(major radius, minor radius)
// p.xyz : rotating the circle around the y-axis
// p.yzx : rotating the circle around the z-axis
// p.zxy : rotating the circle around the x-axis
{
	vec2 q = vec2(length(p.xz) - t.x, p.y);	// zaxis = major axis of rotation
	return length(q) - t.y;
}
float length2( vec2 p )
{
	return sqrt( p.x*p.x + p.y*p.y );
}
float length8( vec2 p )
{
	p = p*p; p = p*p; p = p*p;
	return pow( p.x + p.y, 1.0/8.0 );
}
float fTorus82( vec3 p, vec2 t )	// similar to fTorus
{
	vec2 q = vec2(length2(p.xz) - t.x, p.y);
	return length8(q) - t.y;
}
float fTorus88( vec3 p, vec2 t )	// similar to fTorus
{
	vec2 q = vec2(length8(p.xz) - t.x, p.y);
	return length8(q) - t.y;
}
float fCylinder( vec3 p, vec2 h )
// p = pos relative to center point
// h.x = radius, h.y = half of height
{
	vec2 d = abs( vec2( length(p.xz), p.y ) ) - h;
	return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
}
float fCone( in vec3 p, in vec3 c )
// p = pos relative to apex
// c.xy = normal to cone surface (cf: if c.x>c.y, more steeper cone)
// c.z = height downward from apex
{
	vec2 q = vec2( length(p.xz), p.y );
	float d1 = -q.y-c.z;
	float d2 = max( dot(q,c.xy), q.y );
	return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.0);
}
float fConeSection( in vec3 p, in float h, in float r1, in float r2 )	// h=height, r1=base_radius, r2=top_radius
// p = pos relative to base center
// h = half of height
// r1 = base radius, r2 = top radius
{
	float d1 = -p.y - h;
	float q = p.y - h;
	float si = 0.5*(r1-r2)/h;
	float d2 = max( sqrt( dot(p.xz,p.xz)*(1.0-si*si)) + q*si - r2, q );
	return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.0);
}
float fTriPrism( vec3 p, vec2 h )
// rect-face is parallel to ground & tri-face is normal to ground
// h.x = side length of triangle
// h.y = height of prism (along z-axis)
{
	vec3 q = abs(p);
	return max(q.z-h.y,max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5);
}
float fHexPrism( vec3 p, vec2 h )
// h.x = radius of the inscribed circle of hex
// h.y = half of height (along z-axis)
{
	vec3 q = abs(p);
#if 0
	return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
#else
	float d1 = q.z-h.y;
	float d2 = max((q.x*0.866025+q.y*0.5),q.y)-h.x;
	return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.0);
#endif
}

float fCapsule( vec3 p, vec3 a, vec3 b, float r )
// a = start point, b = end point
// r = radius of both ends
{
	vec3 pa = p - a, ba = b - a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return length( pa - ba*h ) - r;
}
float fPlane( vec3 p, vec3 n )
// p = pos relative to (0,0,0)
// n = unit normal
{
	return dot(n, p);
}

//----------------------------------------------------------------------------
// (function-level) operations (real-valued function: f=(dist to surface))
//----------------------------------------------------------------------------
float fUnion( const in float f1, const in float f2 )
{
	return mix(f1, f2, step(f2, f1));
}
float fIntersect( const in float f1, const in float f2 )
{
	return mix(f2, f1, step(f2, f1));
}
float fSubtract( const in float f1, const in float f2 )
{
	return fIntersect(f1, -f2);
}
float fDisplace( const in float f1, const in float ds )
{
	return f1+ds;
}
float fSmoothPoly( const in float a, const in float b, const in float k ) // polynomial smooth min (k = 0.1)
{
	float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
	return mix( b, a, h ) - k*h*(1.0-h);
}
float fSmoothExp( const in float a, const in float b, const in float k ) // exponential smooth min (k = 32)
{
	float res = exp( -k*a ) + exp( -k*b );
	return -log( res )/k;
}
float fSmoothPow( const in float a, const in float b, const in float k ) // power smooth min (k = 8)
{
	float aa = pow( a, k );
	float bb = pow( b, k );
	return pow( (aa*bb) / (aa+bb), 1.0/k );
}
float fNoise( const in vec2 p )
{
	vec2 s = sin(p * 0.6345) + sin(p * 1.62423);
	return dot(s, vec2(0.125)) + 0.5;
}

//----------------------------------------------------------------------------
// (object-level) operations (vector-valued object: d.xyzw=(dist, matID, texU, texV))
//----------------------------------------------------------------------------
vec4 dUnion( const in vec4 d1, const in vec4 d2 )
{
	return mix(d1, d2, step(d2.x, d1.x));
}
vec4 dUnionTransp( const in vec4 d1, const in vec4 d2, const in float fTranspScale )
{
	// Negate the distance to the transparency object if transparent scale is 0.0
	// This allows us to retrace out of transparency
	vec4 vScaled = vec4(d2.x * (fTranspScale * 2.0 - 1.0), d2.yzw);

	// The condition allows us to ignore transparency for secondary rays
	return mix(d1, vScaled, step(vScaled.x, d1.x) * step(0.0, fTranspScale));
}
vec4 dIntersect( const in vec4 d1, const in vec4 d2 )
{
	return mix(d2, d1, step(d2.x,d1.x));
}
vec4 dSubtract( const in vec4 d1, const in vec4 d2 )
{
	return dIntersect(d1, vec4(-d2.x, d2.yzw));
}
vec4 dDisplace( const in vec4 d1, const in float ds )
// ds = displacement to apply to d1
// for example:
// float ds(vec3 p) { return sin(2*p.x)*sin(2*p.y)*sin(2*p.z); }
// vec4 new_sphere = dDisplace( vec4(fSphere(p,r), gold), ds(p) );
{
	return vec4( d1.x+ds, d1.yzw );
}
vec4 dBlend( const in vec4 d1, const in vec4 d2, const in float k )
// using polynomial smooth min (often k = 0.1)
{
	float fdist = fSmoothPoly( d1.x, d2.x, k );
	vec3 vmat = mix( d1.yzw, d2.yzw, step(d2.x, d1.x) );
	return vec4( fdist, vmat );
}

//----------------------------------------------------------------------------
// (space-level) operations (vector-valued space: p=(x,y,z))
//----------------------------------------------------------------------------
vec3 sRepeat( const in vec3 p, const in vec3 spacing )
// The same domain space is repeated along (x,y,z)-axis, whose size is spacing.xyz, respectively.
// if spacing.y=0, then infinite column is repeated along x-axis & z-axis
{
	//return mod(p, spacing) - 0.5*spacing;
	vec3 q = p - 0.5*spacing;
	return mod(q, spacing) - 0.5*spacing;
}
vec3 sTwist( const in vec3 p, const in float angle )
// twist by angle per unit-length along y-axis with the right-hand rule
{
	float c = cos( angle*p.y );
	float s = sin( angle*p.y );
	mat2 m = mat2( c, -s, s, c );
	vec2 twist = m*p.zx;
	return vec3( twist.y, p.y, twist.x );
}
vec3 sTwistY( const in vec3 p, const in float angle )	// twist along y-axis
{
	return sTwist( p, angle );
}
vec3 sTwistZ( const in vec3 p, const in float angle )	// twist along z-axis
{
	vec3 q = p.yzx;
	q = sTwist( q, angle );
	return q.zxy;
}
vec3 sTwistX( const in vec3 p, const in float angle )	// twist along x-axis
{
	vec3 q = p.zxy;
	q = sTwist( q, angle );
	return q.yzx;
}
vec3 sTranslate( const in vec3 p, const in vec3 t )
{
	return p - t;
}
vec3 sRotateX( const in vec3 p, const in float angle )
{
	float s = sin(angle);
	float c = cos(angle);
	return vec3( p.x, c*p.y+s*p.z, -s*p.y+c*p.z );
}
vec3 sRotateY( const in vec3 p, const in float angle )
{
	float s = sin(angle);
	float c = cos(angle);
	return vec3( c*p.x+s*p.z, p.y, -s*p.x+c*p.z );
}
vec3 sRotateZ( const in vec3 p, const in float angle )
{
	float s = sin(angle);
	float c = cos(angle);
	return vec3( c*p.x+s*p.y, -s*p.x+c*p.y, p.z );
}
// Scale
// e.g., float fScaled = fSphere( p/s, r ) * s
// only possible to scale uniformly
vec3 sCheapBendY( const in vec3 p, float angle )
// bending moment: yaxis, convex up = zaxis, fixed: (0,0,0)
{
	float c = cos( angle*p.x );
	float s = sin( angle*p.x );
	mat2 m = mat2( c, -s, s, c );
	vec2 b = m*p.zx;
	return vec3( b.y, p.y, b.x );
}
vec3 sCheapBendZ( const in vec3 p, float angle )
// bending moment: zaxis, convex up = xaxis, fixed: (0,0,0)
{
	float c = cos( angle*p.y );
	float s = sin( angle*p.y );
	mat2 m = mat2( c, -s, s, c );
	vec2 b = m*p.xy;
	return vec3( b.x, b.y, p.z );
}
vec3 sCheapBendX( const in vec3 p, float angle )
// bending moment: xaxis, convex up: yaxis, fixed: (0,0,0)
{
	float c = cos( angle*p.z );
	float s = sin( angle*p.z );
	mat2 m = mat2( c, -s, s, c );
	vec2 b = m*p.yz;
	return vec3( p.x, b.x, b.y );
}

//----------------------------------------------------------------------------
// Scene Description
//----------------------------------------------------------------------------
vec4 dCheckerBoard( const in vec3 vPos )
{
	return vec4( fPlane( vPos-vec3(0.0), vec3(0.0, 1.0, 0.0) ), kMaterialGround, vPos.xy );
}
vec4 SimpleScene( const in vec3 vPos, const in float fTranspScale )
{
	vec3 tex0 = vec3(kMaterialTexture0, vPos.z, vPos.x);
	vec3 tex1 = vec3(kMaterialTexture1, vPos.x, vPos.y);
	vec3 gold = vec3(kMaterialGold, vPos.z, vPos.x);
	vec3 silver = vec3(kMaterialSilver, vPos.z, vPos.x);
	vec4 res = dCheckerBoard( vPos );
#if 0
	res = dUnion( res, vec4( fBox(vPos-vec3(0.5,0.5,0.5), vec3(0.5)), silver ) );
#endif

#if 0
	res = dUnion( res, vec4( fBox(vPos-vec3(0.5,0.5,0.5), vec3(0.25)), tex0 ) );
	res = dUnion( res, vec4( fCylinder(vPos-vec3(1,0.3,0), vec2(0.2, 0.3)), tex1 ) );
	res = dUnion( res, vec4( fCone(vPos-vec3( 0.0,0.9, 1.0), vec3(0.8, 0.1, 0.6) ), gold ) );
#endif

#if 1
	vec3 p = sRepeat( vPos, vec3(10.0, 10.0, 10.0) );
	res = dUnion( res, vec4( fSphere(p, 2.5), tex0 ) );
#endif

#if 0
	vec3 p = sCheapBendX( vPos, kPI/4.0 );
	res = dUnion( res, vec4( fBox(p, vec3(1, 0.2, 1)), gold ) );
#endif

#if 0
	float a = 2.0;	float r = 0.5;
	vec3 p = sTranslate( vPos, vec3(0,0,0) );
	float dfmValue = sin(a*p.x)*sin(a*p.y)*sin(a*p.z);
	vec4 dsphere = vec4( fSphere(p, r), gold );
	vec4 dsphere1 = dDisplace( dsphere, dfmValue );
	res = dUnion( res, dsphere1 );	//psk : some portion is clipped out (need to debug)
#endif

#if 0
	float k = 0.1;
	vec3 p = sTranslate( vPos, vec3(0,2,0) );
	vec4 dcy1 = vec4( fCylinder(p, vec2(1.0,2.0)), gold );
	vec4 dcy2 = vec4( fCylinder(p.yzx, vec2(0.5, 2.0)), silver );
	res = dUnion( res, dBlend( dcy1, dcy2, k ) );
#endif

	return res;
}
vec4 PrimitiveScene( const in vec3 pos, const in float fTranspScale )
{
	//vec3 gold = vec3(kMaterialGold, 0, 0);
	vec3 gold = vec3(kMaterialGold, pos.x, pos.y);	//psk : texUV=(pos.x, pos.y)
	
	vec4 res = dCheckerBoard( pos );
	res = dUnion( res, vec4( fSphere(    pos-vec3( 0.0,0.25, 0.0), 0.25 ), gold ) );
	res = dUnion( res, vec4( fBox(       pos-vec3( 1.0,0.25, 0.0), vec3(0.25) ), gold ) );
	res = dUnion( res, vec4( fTorus(     pos-vec3( 0.0,0.25, 1.0), vec2(0.20,0.05) ), gold ) );
	res = dUnion( res, vec4( fCapsule(   pos,vec3(-1.3,0.10,-0.1), vec3(-0.8,0.50,0.2), 0.1  ), gold ) );
	res = dUnion( res, vec4( fTriPrism(  pos-vec3(-1.0,0.25,-1.0), vec2(0.25,0.05) ), gold ) );
	res = dUnion( res, vec4( fCylinder(  pos-vec3( 1.0,0.30,-1.0), vec2(0.1,0.2) ), gold ) );
	res = dUnion( res, vec4( fCone(      pos-vec3( 0.0,0.50,-1.0), vec3(0.8,0.6,0.3) ), gold ) );
	res = dUnion( res, vec4( fTorus82(   pos-vec3( 0.0,0.25, 2.0), vec2(0.20,0.05) ), gold ) );
	res = dUnion( res, vec4( fTorus88(   pos-vec3(-1.0,0.25, 2.0), vec2(0.20,0.05) ), gold ) );
	res = dUnion( res, vec4( fCylinder(  pos-vec3( 1.0,0.30, 2.0), vec2(0.1,0.2) ), gold ) );
	res = dUnion( res, vec4( fHexPrism(  pos-vec3(-1.0,0.20, 1.0), vec2(0.25,0.05) ), gold ) );
	res = dUnion( res, vec4( fSubtract(  fBox( pos-vec3(-2.0,0.2, 1.0), vec3(0.2)), fSphere( pos-vec3(-2.0,0.2, 1.0), 0.25)), gold ) );
	res = dUnion( res, vec4( fSubtract( fTorus82(  pos-vec3(-2.0,0.2, 0.0), vec2(0.20,0.1)), fCylinder( sRepeat( vec3(atan(pos.x+2.0,pos.z)/6.2831, pos.y, 0.02+0.5*length(pos-vec3(-2.0,0.2,0.0)) ), vec3(0.05,1.0,0.05)), vec2(0.02,0.6))), gold ) );
	res = dUnion( res, vec4( 0.7*fSphere( pos-vec3(-2.0,0.25,-1.0), 0.2 ) + 0.03*sin(50.0*pos.x)*sin(50.0*pos.y)*sin(50.0*pos.z), gold ) );
	res = dUnion( res, vec4( 0.5*fTorus( sTwist(pos-vec3(-2.0,0.25, 2.0), 20.0), vec2(0.20,0.05)), gold ) );
	res = dUnion( res, vec4( fConeSection( pos-vec3( 0.0,0.35,-2.0), 0.15, 0.2, 0.1 ), gold ) );
	res = dUnion( res, vec4( fEllipsoid( pos-vec3( 1.0,0.35,-2.0), vec3(0.15, 0.2, 0.05) ), gold ) );
	return res;
}
vec4 SpoutScene( const in vec3 vPos, const in float fTranspScale )
{
    const float kPipeRadius = 0.4;
    const float kPipeThickness = 0.15;
    const float kWaterNoiseScale = 0.025;
    const float kWaterVelocity = 1.0;
    const float kWaterAccel = -1.0;
    const float kWaterAnimSpeed = 80.0;
    const float kTrenchWaterAnimSpeed = 20.0;

    //float kPipeHeight = 2.0;
	float kPipeHeight = 2.0 + sin(time);
    float kRipplePos = sqrt(abs(2.0 * kPipeHeight / kWaterAccel)) * kWaterVelocity;

	vec4 vResult = vec4(10000.0, -1.0, 0.0, 0.0);
	float fDistFloor = vPos.y;
	float fDistBrick = fDistFloor;
	float fDistTrench = length(vPos.yz + vec2(-0.4, 0.0)) - 1.0;
	fDistBrick = max(fDistBrick, -(fDistTrench));
	float fDistWall = vPos.x + 1.0;
	fDistBrick = min(fDistBrick, fDistWall);
	vec4 vDistFloor = vec4(fDistBrick, kMaterialWall, vPos.xz + vec2(vPos.y, 0.0));
	vResult = dUnion(vResult, vDistFloor);
	vec3 vWaterDomain = vPos - vec3(0.0, kPipeHeight, 0.0);

	float t = max(vWaterDomain.x / kWaterVelocity, 0.0);

	// Equations of motion
	float s = 0.5 * kWaterAccel * t * t;
	float v = -kWaterAccel * t;

	vWaterDomain.y -= s;
	float fDistWater = (length(vWaterDomain.yz) - kPipeRadius);
	float fDistPipe = max(fDistWater - kPipeThickness, vWaterDomain.x);
	fDistPipe = max(fDistPipe, -fDistWater); // subtract the water from the pipe to make the hole
	vec4 vDistPipe = vec4(fDistPipe, kMaterialPipe, vPos.xy);

	vResult = dUnion(vResult, vDistPipe);

	// compensate for domain distortion of water, otherwise ray sometimes misses
	fDistWater /= (1.0 + v * 0.5);
	vec2 vNoiseDomain = vPos.xz;

	// modify noise for water in trench
	float fInTrench = step(vPos.y, (-0.1 + 0.05));
	vec2 vRippleCentre1 = vPos.xz - vec2(kRipplePos, 0.0);
	vNoiseDomain.x = mix(vNoiseDomain.x, length(vRippleCentre1), fInTrench);
	float fNoiseScale = mix(t * t, 1.0 / (1.0 + vNoiseDomain.x), fInTrench) * kWaterNoiseScale;
	float fWaterSpeed = mix(kWaterAnimSpeed * kWaterVelocity, kTrenchWaterAnimSpeed, fInTrench);

	vNoiseDomain *= 30.0;
	vNoiseDomain.x += -time * fWaterSpeed;

	float fTrenchWaterDist = vPos.y + 0.1;
	fDistWater = min(fDistWater, fTrenchWaterDist);
	fDistWater += fNoise(vNoiseDomain) * fNoiseScale;
	vec4 vDistWater = vec4(fDistWater, kMaterialWater, vPos.xy);
	vResult = dUnionTransp(vResult, vDistWater, fTranspScale);
	return vResult;
}

//----------------------------------------------------------------------------

vec4 GetDistanceScene( const in vec3 vPos, const in float fTranspScale )
// return vec4.x = scene_distance
//        vec4.y = material (or object) id
//        vec4.zw = material specific parameters (maybe uv coordinates)
{
	return SpoutScene( vPos, fTranspScale );		// scene example 1
	//return SimpleScene( vPos, fTranspScale );
	//return PrimitiveScene( vPos, fTranspScale );
}

CMaterial GetObjectMaterial( const in CHitInfo hitInfo )
{
#if 0		// reflectance of a material
mat.cAlbedo = vec3(0.93, 0.88, 0.38);	// Gold
mat.cAlbedo = vec3(0.26, 0.28, 0.26);	// Iridium
mat.cAlbedo = vec3(0.44, 0.435, 0.43);	// Iron
mat.cAlbedo = vec3(0.50, 0.47, 0.36);	// Nickel
mat.cAlbedo = vec3(0.93, 0.80, 0.46);	// Copper
mat.cAlbedo = vec3(0.63, 0.62, 0.57);	// Platinum
mat.cAlbedo = vec3(0.97, 0.97, 0.96);	// Silver
#endif
	CMaterial mat;
	if(hitInfo.vObjectID.x == kMaterialTexture0)		// textureMaps[0]
	{
		mat.cAlbedo = texture( textureMaps[0], hitInfo.vObjectID.yz ).rgb;
		mat.fR0 = (1.0003-0.27049)/(1.0003+0.27049);
		mat.fR0 *= mat.fR0;
		mat.fSmoothness = 1.0;
		mat.fTransparency = 0.0;
	}
	else if(hitInfo.vObjectID.x == kMaterialTexture1)		// textureMaps[1]
	{
		mat.cAlbedo = texture( textureMaps[1], hitInfo.vObjectID.yz ).rgb;
		mat.fR0 = 0.2;//2.0;
		mat.fR0 *= mat.fR0;
		mat.fSmoothness = 1.0;
		mat.fTransparency = 0.0;
	}
	else if(hitInfo.vObjectID.x == kMaterialGround)
	{
		float ncell = 2.0;//5.0;	// ncell = no. cell per unit length (e.g. if ncell=2, then cellsize=0.5)
		float f = mod( floor(ncell*hitInfo.vPos.z) + floor(ncell*hitInfo.vPos.x), 2.0 );
		mat.cAlbedo = 0.4 + 0.1*f*vec3(1.0);
		mat.fR0 = 4.0;//0.8;
		mat.fSmoothness = 1.0;
		mat.fTransparency = 0.0;
	}
	else if(hitInfo.vObjectID.x == kMaterialGold)
	{
		mat.cAlbedo = vec3(1.0, 0.84, 0.0);
		mat.fR0 = (1.0003-0.27049)/(1.0003+0.27049);
		mat.fR0 *= mat.fR0;
		mat.fSmoothness = 1.0;
		mat.fTransparency = 0.0;
	}
	else if(hitInfo.vObjectID.x == kMaterialSilver)
	{
		mat.cAlbedo = vec3(0.75, 0.75, 0.75);
		mat.fR0 = (1.0003-0.15016)/(1.0003+0.15016);
		mat.fR0 *= mat.fR0;
		mat.fSmoothness = 1.0;
		mat.fTransparency = 0.0;
	}
	else if(hitInfo.vObjectID.x == kMaterialWall)
	{
		// floor
		mat.fR0 = 0.02;
	#ifdef ENABLE_TEXTURE
		//hitInfo.vObjectID.yz = fragTexCoord;
		vec3 cTexture = texture(textureMaps[0], hitInfo.vObjectID.yz * 0.25).rgb;

		mat.cAlbedo = cTexture * cTexture;
		mat.fSmoothness = mat.cAlbedo.r;
		mat.fTransparency = 0.0;
	#else
		// Textureless version
		vec2 vTile = step(vec2(0.15), fract(hitInfo.vObjectID.yz));
		float fTile = vTile.x * vTile.y;
		mat.cAlbedo = vec3(1.0) * (fTile * 0.8 + 0.2);
		mat.fSmoothness = 1.0;
	#endif
	}
	else if(hitInfo.vObjectID.x == kMaterialPipe)
	{
		mat.fR0 = 0.8;
		mat.fSmoothness = 1.0;
		mat.cAlbedo = vec3(0.5);
		mat.fTransparency = 0.0;
	}
	else if(hitInfo.vObjectID.x == kMaterialWater)
	{
		mat.fR0 = 0.01;
		mat.fSmoothness = 1.0;
		mat.fTransparency = 1.0;
		mat.fRefractIndex = 1.0 / 1.3330;	// n1(air) / n2(water)
		const float fExtinctionScale = 2.0;
		const vec3 vExtinction = vec3(0.3, 0.7, 0.9);
		mat.cAlbedo = (vec3(1.0) - vExtinction) * fExtinctionScale; // becomes extinction for transparency
	}
	else
	{
		// red color not to define a material
		mat.cAlbedo = vec3(1,0,0);
		mat.fSmoothness = 1.0;
		mat.fTransparency = 0.0;
	}
	return mat;
}

vec3 GetSkyGradient( const in vec3 vDir )
{
	const vec3 cColourTop = vec3(0.7, 0.8, 1.0);    // cSky
	const vec3 cColourHorizon = cColourTop * 0.5;   // cGround
	float fBlend = clamp(vDir.y, 0.0, 1.0);         // w = clamp(vDir.y, 0.0, 1.0)
	return mix(cColourHorizon, cColourTop, fBlend); // skyGradient = cGround*(1-w) + cSky * w
}

CPointLight GetPointLight()
{
	CPointLight result;
	result.vPos = vec3(0.5, 1.0, -2.0);
    //result.cCol = vec3(32.0, 6.0, 1.0) * 10.0;
	result.cCol = vec3(32.0, 6.0, 1.0) * 2.0;
	return result;
}

CDirectionalLight GetDirectionalLight()
{
	CDirectionalLight result;
	result.vDir = normalize(vec3(-0.2, -0.3, 0.5));
    //result.cCol = vec3(8.0, 7.5, 7.0);
	result.cCol = vec3(2.0, 1.5, 1.0);
	return result;
}

vec3 GetAmbientLight(const in vec3 vNormal)
{
	return GetSkyGradient(vNormal); // Lc = skyGradient
}

//----------------------------------------------------------------------------
// Raymarching
//----------------------------------------------------------------------------
vec3 GetSceneNormal( const in vec3 vPos, const in float fTranspScale )
{
	vec3 eps = vec3( kRaymarchEpsilon, 0.0, 0.0 );
	vec3 nor = vec3(
		GetDistanceScene(vPos+eps.xyy, fTranspScale).x - GetDistanceScene(vPos-eps.xyy, fTranspScale).x,
		GetDistanceScene(vPos+eps.yxy, fTranspScale).x - GetDistanceScene(vPos-eps.yxy, fTranspScale).x,
		GetDistanceScene(vPos+eps.yyx, fTranspScale).x - GetDistanceScene(vPos-eps.yyx, fTranspScale).x );
	return normalize(nor);
}
void Raymarch( const in CRay ray, out CHitInfo result, const int maxIter, const float fTranspScale )
// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
{
    result.fDist = ray.fStartDist;
    result.vObjectID.x = 0.0;
    for(int i=0;i<=kRaymarchMaxIter; i++)
    {
        result.vPos = ray.vOrigin + ray.vDir * result.fDist;
        vec4 vSceneDist = GetDistanceScene( result.vPos, fTranspScale );
        result.vObjectID = vSceneDist.yzw;
        // abs allows backward stepping - should only be necessary for non uniform distance functions
        if( (abs(vSceneDist.x) <= kRaymarchEpsilon) || (result.fDist >= ray.fLength) || (i > maxIter) ) break;
        result.fDist = result.fDist + vSceneDist.x;
    }
    if(result.fDist >= ray.fLength)
    {
        result.fDist = kFarClip;
        result.vPos = ray.vOrigin + ray.vDir * result.fDist;
        result.vObjectID.x = 0.0; // maybe sky
    }
}

float GetShadow( const in vec3 vPos, const in vec3 vNormal, const in vec3 vLightDir, const in float fLightDistance )
{
#if defined(ENABLE_HARD_SHADOWS)		// #ifdef ENABLE_HARD_SHADOWS
    CRay shadowRay;
    shadowRay.vDir = vLightDir;
    shadowRay.vOrigin = vPos;
    const float fShadowBias = 0.05;
    shadowRay.fStartDist = fShadowBias / abs(dot(vLightDir, vNormal));
    shadowRay.fLength = fLightDistance - shadowRay.fStartDist;
    CHitInfo shadowIntersect;
    Raymarch(shadowRay, shadowIntersect, 32, kTranspNo);
    float fShadow = step(0.0, shadowIntersect.fDist) * step(fLightDistance, shadowIntersect.fDist );
    return fShadow;
#elif defined(ENABLE_SOFT_SHADOWS)
    const float fShadowBias = 0.05;
    float fDist = fShadowBias / abs(dot(vLightDir, vNormal));
    float fDistMax = fLightDistance - fDist;
    float fResult = 1.0;
    for( int i=0; i<16; i++ )
    {
        float fSceneDist = GetDistanceScene( vPos + vLightDir*fDist, kTranspNo ).x;
        fResult = min( fResult, 8.0*fSceneDist/fDist );
        fDist += clamp( fSceneDist, 0.02, 0.10 );			//psk
        if( fSceneDist<0.001 || fDist>fDistMax ) break;
    }
    return clamp( fResult, 0.0, 1.0 );
#else
    return 1.0;
#endif
}

float GetAmbientOcclusion(const in CHitInfo intersection, const in CSurface surface)
{
#ifdef ENABLE_AMBIENT_OCCLUSION
    float fScaleOcc = 0.2;		//psk
    vec3 vPos = intersection.vPos;
    vec3 vNormal = surface.vNormal;
    float fAmbientOcclusion = 1.0;
    float fDist = 0.0;
    for(int i=0; i<=5; i++)
    {
        fDist += 0.1;
        vec4 vSceneDist = GetDistanceScene(vPos + vNormal * fDist, kTranspNo);
        fAmbientOcclusion *= 1.0 - max(0.0, (fDist - vSceneDist.x) * fScaleOcc / fDist );
    }
    return fAmbientOcclusion;
#else
    return 1.0;
#endif
}

//----------------------------------------------------------------------------
// Lighting and Shading
//----------------------------------------------------------------------------
void ApplyAtmosphere(inout vec3 col, const in CRay ray, const in CHitInfo hitInfo)
{
#ifdef ENABLE_FOG
    float fFogAmount = exp(hitInfo.fDist * -kFogDensity);
    vec3 cFog = GetSkyGradient(ray.vDir);
    #ifdef ENABLE_DIRECTIONAL_LIGHT_FLARE
        CDirectionalLight directionalLight = GetDirectionalLight();
        float fDirDot = clamp(dot(-directionalLight.vDir, ray.vDir), 0.0, 1.0);
        cFog += directionalLight.cCol * pow(fDirDot, 10.0);
    #endif
    col = mix(cFog, col, fFogAmount);
#endif
#ifdef ENABLE_POINT_LIGHT_FLARE
    CPointLight pointLight = GetPointLight();
    vec3 vToLight = pointLight.vPos - ray.vOrigin;
    float fPointDot = dot(vToLight, ray.vDir);
    fPointDot = clamp(fPointDot, 0.0, hitInfo.fDist);
    vec3 vClosestPoint = ray.vOrigin + ray.vDir * fPointDot;
    float fDist = length(vClosestPoint - pointLight.vPos);
    col += pointLight.cCol * 0.01/ (fDist * fDist);
#endif
}

float Schlick( const in vec3 vHalf, const in vec3 vView, const in float fR0, const in float fSmoothFactor)
{
	float fDot = dot(vHalf, -vView); // H o V
	fDot = clamp((1.0 - fDot), 0.0, 1.0); // 1 - H o V
	float fDotPow = pow(fDot, 5.0); // (1 - H o V)^5
	return fR0 + (1.0 - fR0) * fDotPow * fSmoothFactor; // fresnel = fR0 + (1.0 - fR0) * (1 - H o V)^5 * smoothFactor
}

vec3 ApplyFresnel(const in vec3 vDiffuse, const in vec3 vSpecular, const in vec3 vNormal, const in vec3 vView, const in CMaterial material)
// vView = ray.vDir
{
	vec3 vReflect = reflect(vView, vNormal); // R = reflect(-V, N)
	vec3 vHalf = normalize(vReflect + -vView); // H = R + V
	float fFresnel = Schlick(vHalf, vView, material.fR0, material.fSmoothness * 0.9 + 0.1);
	return mix(vDiffuse, vSpecular, fFresnel); // return (Ld*Md) * (1-fresnel) + Ls * (fresnel)
}

float GetBlinnPhongIntensity(const in vec3 vIncidentDir, const in vec3 vLightDir, const in vec3 vNormal, const in float fSmoothness)
// vIncidentDir = ray.vDir = eyePos to surfPos
// vLightDir = surfPos to lightPos
{          
	vec3 vHalf = normalize(vLightDir - vIncidentDir);  // H = L + V
	float fNdotH = max(0.0, dot(vHalf, vNormal));      // NoH = N o H
	float fSpecPower = exp2(4.0 + 6.0 * fSmoothness);  // fSmoothness ==> fSpecPower (= shininess)
	float fSpecIntensity = (fSpecPower + 2.0) * 0.125; // fSpecPower ==> fSpecIntensity (= specIntensity)
	return pow(fNdotH, fSpecPower) * fSpecIntensity;   // specWt = specIntensity * NoH^(shininess)
}

CShading ApplyPointLight( const in CPointLight light, const in vec3 vSurfacePos, const in vec3 vIncidentDir, const in vec3 vNormal, const in CMaterial material )
// vIncidentDir = ray.vDir = eyePos to surfPos
{
	CShading shading;
	vec3 vToLight = light.vPos - vSurfacePos; // L
	vec3 vLightDir = normalize(vToLight); // L = L/|L|
	float fLightDistance = length(vToLight); // d = |L|
	float fAttenuation = 1.0 / (fLightDistance * fLightDistance); // atten = 1/(d*d)
	float fShadowFactor = GetShadow( vSurfacePos, vNormal, vLightDir, fLightDistance ); // shadow
	vec3 vIncidentLight = light.cCol * fShadowFactor * fAttenuation * max(0.0, dot(vLightDir, vNormal)); // Ld = shadow * atten * Lc * NoL
	shading.cDiffuse = vIncidentLight; // Ld
	shading.cSpecular = GetBlinnPhongIntensity( vIncidentDir, vLightDir, vNormal, material.fSmoothness ) * vIncidentLight; // Ls = Ld * specWt
	return shading;
}

CShading ApplyDirectionalLight( const in CDirectionalLight light, const in vec3 vSurfacePos, const in vec3 vIncidentDir, const in vec3 vNormal, const in CMaterial material )
// vIncidentDir = ray.vDir = eyePos to surfPos
{
	CShading shading;
	const float kShadowRayLength = 10.0;
	vec3 vLightDir = -light.vDir; // L
	float fShadowFactor = GetShadow( vSurfacePos, vNormal, vLightDir, kShadowRayLength ); // shadow
	vec3 vIncidentLight = light.cCol * fShadowFactor * max(0.0, dot(vLightDir, vNormal)); // Ld = shadow * Lc * NoL
	shading.cDiffuse = vIncidentLight; // Ld
	shading.cSpecular = GetBlinnPhongIntensity( vIncidentDir, vLightDir, vNormal, material.fSmoothness ) * vIncidentLight; // Ls = Ld * specWt
	return shading;
}

vec3 ShadeSurface(const in CRay ray, const in CHitInfo hitInfo, const in CSurface surface, const in CMaterial material)
{
	vec3 cScene;
	CShading shading;
	shading.cDiffuse = vec3(0.0);  // totalLd = 0.0
	shading.cSpecular = vec3(0.0); // totalLs = 0.0
	float fAmbientOcclusion = GetAmbientOcclusion(hitInfo, surface);
	vec3 vAmbientLight = GetAmbientLight(surface.vNormal) * fAmbientOcclusion; // La = Lc(= skyGradient) * AO
	shading.cDiffuse += vAmbientLight;        // totalLd += La
	shading.cSpecular += surface.cReflection; // totalLs += (Ms)

#ifdef ENABLE_POINT_LIGHT
    CPointLight pointLight = GetPointLight(); 
    CShading pointLighting = ApplyPointLight(pointLight, hitInfo.vPos, ray.vDir, surface.vNormal, material);
    shading.cDiffuse += pointLighting.cDiffuse;   // totalLd += Ld
    shading.cSpecular += pointLighting.cSpecular; // totalLs += Ls
#endif

#ifdef ENABLE_DIRECTIONAL_LIGHT
    CDirectionalLight directionalLight = GetDirectionalLight();
    CShading directionLighting = ApplyDirectionalLight(directionalLight, hitInfo.vPos, ray.vDir, surface.vNormal, material);
    shading.cDiffuse += directionLighting.cDiffuse;	  // totalLd += Ld
    shading.cSpecular += directionLighting.cSpecular; // totalLs += Ls
#endif

	vec3 vDiffuseReflection = shading.cDiffuse * material.cAlbedo; // totalLd * (Md)
	vDiffuseReflection = mix(vDiffuseReflection, surface.cTransmission, material.fTransparency);// totalDiff = (totalLd*Md)*(transp-1) + Mt*transp

#ifdef ENABLE_SPECULAR
    cScene = ApplyFresnel(vDiffuseReflection, shading.cSpecular, surface.vNormal, ray.vDir, material);// cScene = (totalDiff) * (1-fresnel) + totalLs * (fresnel)
#else
    cScene = vDiffuseReflection; // cScene = (totalDiff)
#endif

	return cScene;
}

vec3 GetSceneColourSecondary( const in CRay ray );

vec3 GetReflection( const in CRay ray, const in CHitInfo hitInfo, const in CSurface surface )
{
#ifdef ENABLE_REFLECTIONS    
{
	// get colour from reflected ray
    const float fSeparation = 0.1;//psk
    CRay reflectRay;
    reflectRay.vDir = reflect(ray.vDir, surface.vNormal);
    reflectRay.vOrigin = hitInfo.vPos;
    reflectRay.fLength = 16.0;//psk
    reflectRay.fStartDist = fSeparation / abs(dot(reflectRay.vDir, surface.vNormal));
    return GetSceneColourSecondary(reflectRay);
}
#else
    return GetSkyGradient(reflect(ray.vDir, surface.vNormal));
#endif
}

vec3 GetTransmission( const in CRay ray, const in CHitInfo hitInfo, const in CSurface surface, const in CMaterial material )
{
	#ifdef ENABLE_TRANSPARENCY
	{
		const float fSeparation = 0.05;//psk
		// Trace until outside transparent object
		CRay refractRay;
		// we dont handle total internal reflection (in that case refract returns a zero length vector)
		refractRay.vDir = refract(ray.vDir, surface.vNormal, material.fRefractIndex);
		refractRay.vOrigin = hitInfo.vPos;
		refractRay.fLength = 16.0;//psk
		refractRay.fStartDist = fSeparation / abs(dot(refractRay.vDir, surface.vNormal));

		#ifdef DOUBLE_SIDED_TRANSPARENCY
			CHitInfo hitInfo2;
			Raymarch(refractRay, hitInfo2, 32, kTranspInverse);
			vec3 vNormal = GetSceneNormal(hitInfo2.vPos, kTranspInverse);
			
			// get colour from rest of scene
			CRay refractRay2;
			refractRay2.vDir = refract(refractRay.vDir, vNormal, 1.0 / material.fRefractIndex);
			refractRay2.vOrigin = hitInfo2.vPos;
			refractRay2.fLength = 16.0;
			refractRay2.fStartDist = 0.0;//fSeparation / abs(dot(refractRay2.vDir, vNormal));
			
			float fExtinctionDist = hitInfo2.fDist;
			vec3 vSceneColour = GetSceneColourSecondary(refractRay2);
		#else
			vec3 vSceneColour = GetSceneColourSecondary(refractRay);
			float fExtinctionDist = 0.5;//psk
		#endif

		vec3 cMaterialExtinction = material.cAlbedo;
		// extinction should really be exp(-) but this is a nice hack to get RGB
		vec3 cExtinction = (1.0 / (1.0 + (cMaterialExtinction * fExtinctionDist)));
		return vSceneColour * cExtinction;
	}
	#else
		return GetSkyGradient(reflect(ray.vDir, surface.vNormal));
	#endif
}

vec3 GetSceneColourSecondary( const in CRay ray )
// no reflections, no transparency, used for secondary rays
{
	CHitInfo hitInfo;
	Raymarch(ray, hitInfo, 32, kTranspNo);
	vec3 cScene;
	if( hitInfo.vObjectID.x < 0.5 )
	{
		cScene = GetSkyGradient(ray.vDir);
	}
	else
	{
		CSurface surface;
		surface.vNormal = GetSceneNormal(hitInfo.vPos, kTranspNo);
		CMaterial material = GetObjectMaterial(hitInfo);
		// use sky gradient instead of reflection
		surface.cReflection = GetSkyGradient(reflect(ray.vDir, surface.vNormal));
		material.fTransparency = 0.0;
		// apply lighting
		cScene = ShadeSurface(ray, hitInfo, surface, material);
	}
	ApplyAtmosphere(cScene, ray, hitInfo);
	return cScene;
}

vec3 GetSceneColourTestVersion( const in CRay ray )
{
	vec3 cScene;
	CHitInfo intersection;
	Raymarch(ray, intersection, 256, kTranspNo);
	if( intersection.vObjectID.x < 0.5 )
	{
		cScene = GetSkyGradient(ray.vDir);
	}
	else
    {
		float t = intersection.fDist;
        vec3 pos = intersection.vPos;
		CSurface surface;
		surface.vNormal = GetSceneNormal(intersection.vPos, kTranspNo);
		vec3 nor = surface.vNormal;
        vec3 ref = reflect( ray.vDir, nor );
		//vec3 col = 0.45 + 0.3*sin( vec3(0.05,0.08,0.10)*(m-1.0) ); // material
		CMaterial material = GetObjectMaterial( intersection );
		vec3 col = material.cAlbedo;
		float occ = GetAmbientOcclusion( intersection, surface );
		// lighitng
		vec3 lig = normalize( vec3(-0.6, 0.7, -0.5) );
		float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0 );
        float dif = clamp( dot( nor, lig ), 0.0, 1.0 );
        float bac = clamp( dot( nor, normalize(vec3(-lig.x,0.0,-lig.z))), 0.0, 1.0 )*clamp( 1.0-pos.y,0.0,1.0);
        float dom = smoothstep( -0.1, 0.1, ref.y );
        float fre = pow( clamp(1.0+dot(nor, ray.vDir), 0.0,1.0), 2.0 );
		float spe = pow(clamp( dot( ref, lig ), 0.0, 1.0 ),16.0);
		dif *= GetShadow( pos, nor, lig, 2.5 );		// 2.5 = max distance from hitPos to lightPos
		dom *= GetShadow( pos, nor, ref, 2.5 );
		vec3 lin = vec3(0.0);
		float scale = 2.0;//1.2;
        lin += scale*dif*vec3(1.00,0.85,0.55);
		lin += scale*spe*vec3(1.00,0.85,0.55)*dif;
        lin += 0.20*amb*vec3(0.50,0.70,1.00)*occ;
        lin += 0.30*dom*vec3(0.50,0.70,1.00)*occ;
        lin += 0.30*bac*vec3(0.25,0.25,0.25)*occ;
        lin += 0.40*fre*vec3(1.00,1.00,1.00)*occ;
		col = col*lin;
    	col = mix( col, vec3(0.8,0.9,1.0), 1.0-exp( -0.0002*t*t ) ); // fog
		col = clamp(col, 0.0, 1.0);
		cScene = pow( col, vec3(0.4545) ); //psk
    }
	return cScene;
}

vec3 GetSceneColourPrimary( const in CRay ray )
{
#ifdef ENABLE_TEST_RENDER
	return GetSceneColourTestVersion( ray );
#else

	CHitInfo intersection;
	Raymarch(ray, intersection, 256, kTranspYes);
	vec3 cScene;
	if( intersection.vObjectID.x < 0.5 )
	{
		cScene = GetSkyGradient(ray.vDir);

		// fragment depth
		gl_FragDepth = 0.99;
	}
	else
	{
		CSurface surface;
		// surface normal
		surface.vNormal = GetSceneNormal(intersection.vPos, kTranspYes);
		// material selection
		CMaterial material = GetObjectMaterial(intersection);
		// surface reflection
		surface.cReflection = GetReflection(ray, intersection, surface);
		// surface transmission
		if(material.fTransparency > 0.0)
		{
			surface.cTransmission = GetTransmission(ray, intersection, surface, material);
		}
		// apply lighting
		cScene = ShadeSurface(ray, intersection, surface, material);

		// fragment depth
		gl_FragDepth = getFragDepth( intersection.vPos );
	}
	ApplyAtmosphere(cScene, ray, intersection);
	return cScene;
#endif
}

vec3 Tonemap( const in vec3 cCol )
{
	vec3 vResult = 1.0 -exp2(-cCol);
	return vResult;
}

void main()
{
	@import ./ray;

	CRay ray;
	ray.vOrigin = ro;
	ray.vDir = rd;
	ray.fStartDist = 0.0;
	ray.fLength = kFarClip;

	vec3 cScene = GetSceneColourPrimary( ray );

	gl_FragColor = vec4( Tonemap(cScene * 1.5), 1.0 ); //fExposure = 1.5
}