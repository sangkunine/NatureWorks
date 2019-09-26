@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;
@import ./scene;

//==============================================================================
// Enables
//==============================================================================

//#define ENABLE_TEST_RENDER

// Choose...
#define ENABLE_SPECULAR
#define ENABLE_REFLECTIONS
// ENABLE_SPECULAR(o) && ENABLE_REFLECTIONS(o) ==> glossy, mirror
// ENABLE_SPECULAR(o) && ENABLE_REFLECTIONS(x) ==> glossy only
// ENABLE_SPECULAR(x) && ENABLE_REFLECTIONS(o) ==> no glossy, no mirror
// ENABLE_SPECULAR(x) && ENABLE_REFLECTIONS(x) ==> no glossy, no mirror
#define ENABLE_TRANSPARENCY
#define ENABLE_DOUBLE_TRANSPARENCY
#define ENABLE_FOG
#define ENABLE_BUMP_MAPPING
//#define ENABLE_SUBSURFACE		// <=== not complete...

// Choose...
#define ENABLE_SHADOW
#define ENABLE_AO

// Choose...
#define ENABLE_SUNDIR_LIGHT
#define ENABLE_DIRECTIONAL_LIGHT
//#define ENABLE_POINT_LIGHT
//#define ENABLE_DIRECTIONAL_LIGHT_FLARE
//#define ENABLE_POINT_LIGHT_FLARE

// Choose...
//#define DEMO_REPEAT_SPHERES
//#define DEMO_PRIMITIVES
#define DEMO_TEST

#ifdef DEMO_PRIMITIVES
	#undef ENABLE_TRANSPARENCY
	#undef ENABLE_DOUBLE_TRANSPARENCY
	#undef ENABLE_FOG
	#undef ENABLE_REFLECTIONS
	#undef ENABLE_SPECULAR
#endif

//==============================================================================
// constants
//==============================================================================

const vec3 SUN_LIGHT = normalize(vec3(0.4, 0.4, -0.48));
const vec3 SUN_COLOR = vec3(1.0, 0.9, 0.83);

const float FAR = 100.0;
const float airIOR = 1.0003;
const float fogDensity = 0.001;//0.05

const float kMaterialSky = 0.0;
const float kMaterialChecker = 1.0;
const float kMaterialGround = 2.0;
const float kMaterialGold = 3.0;
const float kMaterialSilver = 4.0;
const float kMaterialTile = 5.0;
const float kMaterialPipe = 6.0;
const float kMaterialWater = 7.0;
const float kMaterialTexture0 = 11.0;	// textureMaps[0]
const float kMaterialTexture1 = 12.0;	// textureMaps[1]
const float kMaterialTexture2 = 13.0;	// textureMaps[1]
const float kMaterialTexture3 = 14.0;	// textureMaps[1]

const vec3 cGold = vec3(1.000, 0.766, 0.336);
const vec3 cSilver = vec3(0.972, 0.960, 0.915);
const vec3 cCopper = vec3(0.955, 0.637, 0.538);
const vec3 cIron = vec3(0.56, 0.57, 0.58);
const vec3 cAluminium = vec3(0.913, 0.921, 0.925);
const vec3 cChromium = vec3(0.550, 0.556, 0.554);
const vec3 cNickel = vec3(0.660, 0.609, 0.526);
const vec3 cTitanium = vec3(0.542, 0.497, 0.449);
const vec3 cCobalt = vec3(0.662, 0.655, 0.634);
const vec3 cPlatinum = vec3(0.672, 0.637, 0.585);

//==============================================================================
// Data structures
//==============================================================================

struct CRay
{
	vec3 o;
	vec3 d;
	float tmin;
	float tmax;	// maximum distance of ray
};
struct CHit
{
	vec3 p;		// p = ro + rd * dist
	vec3 n;		// surface normal
	float t;
	float m;
	vec3 refCol; // reflected color
	vec3 trmCol; // transmitted color (refraction)
};
struct CMaterial
{
	vec3 albedo;		// diffuse reflectivity
	float ior;			// F0 = ((1-ior)/(1+ior))^2 (where ior = Index of Refraction)
	float smoothness;	// Fresnel = F0 + (1.0 - F0) * (1 - cosTheta)^5 * smoothness
	float transparency; // diffuse*(1-fTransparency) + cTransmission*fTransparency
};
struct CShading
{
	vec3 diffuse;
	vec3 specular;
};
struct CPointLight
{
	vec3 pos;
	vec3 col;
};
struct CDirLight
{
	vec3 dir;
	vec3 col;
};

//==============================================================================
// Scene Description
//==============================================================================

#ifdef DEMO_REPEAT_SPHERES
vec2 demoRepeatSpheres( in vec3 p )
{
	vec2 d = vec2( fPlane(p), kMaterialChecker ); // base plane

	vec3 p1 = sRepeat( p, vec3(10.0, 10.0, 10.0) );
	d = dUnion( d, vec2( fSphere(p1, 2.5), kMaterialSilver ) );//kMaterialTexture0

	return d;
}
#endif

#ifdef DEMO_PRIMITIVES
vec2 demoPrimitives(vec3 p)
{
    // 21 objects + 1 floor
    float m = 0.0;
    vec2 res = dUnion( vec2( fPlane(     p), kMaterialChecker ),
	                   vec2( fSphere(    p-vec3( 0.0,0.25, 0.0), 0.25 ), m ) ); m += 1.0;
    res = dUnion( res, vec2( fBox(       p-vec3( 1.0,0.25, 0.0), vec3(0.25) ), m ) ); m += 1.0;
    res = dUnion( res, vec2( uRoundBox(  p-vec3( 1.0,0.25, 1.0), vec3(0.15), 0.1 ), m ) ); m += 1.0;
	res = dUnion( res, vec2( fTorus(     p-vec3( 0.0,0.25, 1.0), vec2(0.20,0.05) ), m ) );  m += 1.0;
    res = dUnion( res, vec2( fCapsule(   p,vec3(-1.3,0.10,-0.1), vec3(-0.8,0.50,0.2), 0.1  ), m ) ); m += 1.0;
	res = dUnion( res, vec2( fTriPrism(  p-vec3(-1.0,0.25,-1.0), vec2(0.25,0.05) ), m ) ); m += 1.0;
	res = dUnion( res, vec2( fCylinder(  p-vec3( 1.0,0.30,-1.0), vec2(0.1,0.2) ), m ) ); m += 1.0;
	res = dUnion( res, vec2( fCone(      p-vec3( 0.0,0.50,-1.0), vec3(0.8,0.6,0.3) ), m ) ); m += 1.0;
	res = dUnion( res, vec2( fTorus82(   p-vec3( 0.0,0.25, 2.0), vec2(0.20,0.05) ), m ) ); m += 1.0;
	res = dUnion( res, vec2( fTorus88(   p-vec3(-1.0,0.25, 2.0), vec2(0.20,0.05) ), m ) ); m += 1.0;
	res = dUnion( res, vec2( fCylinder6( p-vec3( 1.0,0.30, 2.0), vec2(0.1,0.2) ), m ) );  m += 1.0;
	res = dUnion( res, vec2( fHexPrism(  p-vec3(-1.0,0.20, 1.0), vec2(0.25,0.05) ), m ) ); m += 1.0;
	res = dUnion( res, vec2( fPryamid4(  p-vec3(-1.0,0.15,-2.0), vec3(0.8,0.6,0.25) ), m ) );  m += 1.0;
    res = dUnion( res, dSubtract( vec2(  uRoundBox(p-vec3(-2.0,0.2, 1.0), vec3(0.15),0.05), m),
	                              vec2(  fSphere(p-vec3(-2.0,0.2, 1.0), 0.25), m+1.0) ) ); m += 2.0;
    res = dUnion( res, dSubtract( vec2(  fTorus82(p-vec3(-2.0,0.2, 0.0), vec2(0.20,0.1)), m),
	                              vec2(  fCylinder(sRepeat(vec3(atan(p.x+2.0,p.z)/6.2831, p.y, 0.02+0.5*length(p-vec3(-2.0,0.2, 0.0))), vec3(0.05,1.0,0.05)), vec2(0.02,0.6)), m+1.0) )); m += 2.0;
	res = dUnion( res, vec2( 0.5*fSphere( p-vec3(-2.0,0.25,-1.0), 0.2 ) + 0.03*sin(50.0*p.x)*sin(50.0*p.y)*sin(50.0*p.z), m ) ); m += 1.0;
	res = dUnion( res, vec2( 0.5*fTorus(  sTwist(p-vec3(-2.0,0.25,2.0), 10.0), vec2(0.20,0.05)), m ) ); m += 1.0;
    res = dUnion( res, vec2( fConeSection(p-vec3( 0.0,0.35,-2.0), 0.15, 0.2, 0.1 ), m ) ); m += 1.0;
    res = dUnion( res, vec2( fEllipsoid(  p-vec3( 1.0,0.35,-2.0), vec3(0.15, 0.2, 0.05) ), m ) ); m += 1.0;
    return res;
}
#endif

#ifdef DEMO_TEST
vec2 demoTest( in vec3 p )
{
	// // case 1: static sphere
	// vec2 d = vec2( fPlane(p), kMaterialChecker ); // base plane
	// d = dUnion( d, vec2( fSphere(p-vec3(0.0, 2.5, 0.0), 2.5), kMaterialWater ) );

	// // case 2: smooth union/subtract
	// vec2 d = vec2( fPlane(p), kMaterialChecker ); // base plane
	// vec3 tp = sTranslate( p, vec3(0.0, 2.0, 0.0) );
	// vec2 dcy1 = vec2( fCylinder(tp, vec2(1.0, 2.0)), kMaterialTexture0);
	// vec2 dcy2 = vec2( fCylinder(tp.yzx, vec2(0.5, 2.0)), kMaterialTexture1 );
	// d = dUnion( d, dSmoothUnion( dcy1, dcy2, 0.1 ) );
    // //d = dUnion( d, dSmoothSubtract( dcy1, dcy2, 0.1 ) );

	// case 3: wriggling ice
	vec2 d = vec2( fPlane(p), kMaterialChecker ); // base plane
	vec2 sph = vec2( fSphere(p-vec3(0.0, 2.5, 0.0), 2.5), kMaterialTexture3 );
	p = rotM3*p*0.12; sph.x += 0.13*sinusoidBumps( p, time*1.1 );
	p = rotM3*p*0.57; sph.x -= 0.19*sinusoidBumps( p, time*2.1 );
	d = dUnion( d, sph );

	// // case 4: heart
	// vec2 d = vec2( fPlane(p), kMaterialChecker ); // base plane
	// vec2 sph = vec2( fSphere(p-vec3(0.0, 2.5, 0.0), 2.5), kMaterialTexture2 );
	// p = rotM3*p*0.15; sph.x += 0.23*sinusoidBumps( p, time*0.1 );
	// p = rotM3*p*0.15; sph.x -= 0.05*sinusoidBumps( p, time*3.1 );
	// d = dUnion( d, sph );

	return d;
}
#endif

vec2 sceneMap( in vec3 p )
{
#if defined (DEMO_REPEAT_SPHERES)
	return demoRepeatSpheres( p );

#elif defined (DEMO_PRIMITIVES)
	return demoPrimitives( p );

#elif defined (DEMO_TEST)
	return demoTest( p );

#endif
}

//==============================================================================
// Materials
//==============================================================================

CMaterial getMaterial( in CHit hit )
{
	// IOR tables
	// Acetone 1.36
	// Actinolite 1.618
	// Agalmatoite 1.550
	// Agate 1.544
	// Agate, Moss 1.540
	// Air 1.0002926 (= 1.0003)
	// Alcohol 1.329
	// Amber 1.546
	// Amethyst 1.544
	// Crystal 2.00
	// Diamond 2.417
	// Emerald 1.576
	// Ethanol 1.36
	// Ethyl Alcohol 1.36
	// Glass 1.51714
	// Glass, Albite 1.4890
	// Glass, Crown 1.520
	// Glass, Crown, Zinc 1.517
	// Glass, Flint, Dense 1.66
	// Glass, Flint, Heaviest 1.89
	// Glass, Flint, Heavy 1.65548
	// Glass, Flint, Lanthanum 1.80
	// Glass, Flint, Light 1.58038
	// Glass, Flint, Medium 1.62725
	// Gold 0.47
	// Ice 1.309
	// Ivory 1.540
	// Jade, Nephrite 1.610
	// Jadeite 1.665
	// Lead 2.01
	// Malachite 1.655
	// Methanol 1.329
	// Moonstone, Albite 1.535
	// Nylon 1.53
	// Onyx 1.486
	// Opal 1.450
	// Oxygen (gas) 1.000276
	// Oxygen (liq) 1.221
	// Pearl 1.530
	// Plastic 1.460
	// Plexiglas 1.50
	// Polystyrene 1.55
	// Quartz 1.544
	// Quartz, Fused 1.45843
	// Rock Salt 1.544
	// Rubber, Natural 1.5191
	// Ruby 1.760
	// Sapphire 1.760
	// Silicon 4.24
	// Steel 2.50
	// Tiger eye 1.544
	// Topaz 1.620
	// Tourmaline 1.624
	// Turpentine 1.472
	// Turquoise 1.610
	// Water (gas) 1.000261
	// Water 35'C (Room temp) 1.33157
	// Zirconia, Cubic 2.170

	// Silver 0.180
	// Iron 2.950
	// Milk 1.350
	// Copper 1.100
	// Bronze1.180
	// Asphalt 1.635
	// Eye,Lens 1.410
	// Leather, Stone(Pebbles), Wood, Plaster 1.52

	CMaterial m;
	if(hit.m == kMaterialTexture0)		// textureMaps[0]
	{
		m.albedo = texCube( textureMaps[0], hit.p, hit.n ).rgb;
		m.ior = 1.52;//Pebbles
		m.smoothness = 0.0;
		m.transparency = 0.0;
	}
	else if(hit.m == kMaterialTexture1)		// textureMaps[1]
	{
		m.albedo = texCube( textureMaps[1], hit.p, hit.n ).rgb;
		m.ior = 2.50;//Steel
		m.smoothness = 1.0;
		m.transparency = 0.0;
	}
	else if(hit.m == kMaterialTexture2)		// textureMaps[2]
	{
		m.albedo = texCube( textureMaps[2], hit.p, hit.n ).rgb;
		m.ior = 1.760;//Ruby
		m.smoothness = 1.0;
		m.transparency = 0.0;
	}
	else if(hit.m == kMaterialTexture3)		// textureMaps[3]
	{
		m.albedo = texCube( textureMaps[3], hit.p, hit.n ).rgb;
		m.ior = 1.309;//Ice
		m.smoothness = 1.0;
		m.transparency = 1.0;
	}
	else if(hit.m == kMaterialChecker)
	{
		float f = checkerGradBox( 1.0*hit.p.xz );
		m.albedo = 0.1 + f*vec3(0.8);
		//m.ior = 1.309;//Ice
		m.ior = 1.0;//No_reflection
		m.smoothness = 1.0;
		m.transparency = 0.0;
	}
	else if(hit.m == kMaterialGround)
	{
		float f = checkerGradBox( 1.0*hit.p.xz );
		m.albedo = cCopper + 0.1*f*vec3(1.0);
		m.ior = 1.0;//No_reflection
		m.smoothness = 0.0;
		m.transparency = 0.0;
	}
	else if(hit.m == kMaterialGold)
	{
		m.albedo = cGold;
		m.ior = 0.47;//Gold
		m.smoothness = 1.0;
		m.transparency = 0.0;
	}
	else if(hit.m == kMaterialSilver)
	{
		m.albedo = cSilver;
		m.ior = 0.180;//Silver
		m.smoothness = 1.0;
		m.transparency = 0.0;
	}
	else if(hit.m == kMaterialTile)
	{
		// Textureless version
		vec2 vTile = step(vec2(0.15), fract(hit.p.xz));
		float fTile = vTile.x * vTile.y;
		m.albedo = vec3(1.0) * (fTile * 0.8 + 0.2);
		m.ior = 2.01;//Lead
		m.smoothness = m.albedo.r;
		m.transparency = 0.0;
	}
	else if(hit.m == kMaterialPipe)
	{
		m.albedo = vec3(0.5);
		m.ior = 2.50;//Steel
		m.smoothness = 1.0;
		m.transparency = 0.0;
	}
	else if(hit.m == kMaterialWater)
	{
		m.ior = 1.3330; // Water at 20'C
		m.smoothness = 1.0;
		m.transparency = 1.0;
		float fExtinctionScale = 2.0;
		//vec3 vExtinction = vec3(0.3,0.7,0.9);
		vec3 vExtinction = vec3(1.0);
		m.albedo = (vec3(1.0) - vExtinction) * fExtinctionScale;
	}
	else
	{
		m.albedo = hsl2rgb( vec3(hit.m/7.0, 1.0, 0.5) );
		m.ior = 1.410;//Eye,Lens
		m.smoothness = 1.0;
		m.transparency = 0.0;
	}
	return m;
}

//==============================================================================
// Lights
//==============================================================================

CPointLight getPointLight()
{
	CPointLight res;
	res.pos = vec3(2.0, 4.0, -2.0);
	res.col = vec3(1.0, 0.0, 0.0);
	return res;
}

CDirLight getDirLight()
{
	CDirLight res;
	res.dir = vec3(1.0);
	res.col = vec3(1.0);
	return res;
}

CDirLight getSunLight()
{
	CDirLight res;
	res.dir = SUN_LIGHT;
	res.col = SUN_COLOR;
	return res;
}

vec3 getSkyGradient( in vec3 rd )
{
#if 1
	float sunAmount = max(dot(rd, SUN_LIGHT), 0.0);
	float v = pow(1.0 - max(rd.y, 0.0), 5.0)*0.5;
	vec3 sky = v * SUN_COLOR*vec3(0.4) + vec3(0.18,0.22,0.4);
	// wide glare effect...
	sky = sky + SUN_COLOR * min(pow(sunAmount, 60.5)*0.32, 0.3);
	sky = sky+ SUN_COLOR * min(pow(sunAmount, 1150.0), 0.3)*0.65;
	return sky;
#else
	return vec3(0.1);
#endif
}

vec3 getAmbientLight( in vec3 n )
{
	return getSkyGradient( n );
}

//==============================================================================
// Raymarching
//==============================================================================

void rayMarching( in CRay r, out CHit hit )
{
	vec2 tm = rayMarching( r.o, r.d, r.tmin, FAR );
	hit.t = tm.x;
	hit.p = r.o + r.d * hit.t;
	hit.m = (tm.y == -1.0)? kMaterialSky : tm.y;
}

//==============================================================================
// Lighting & Shadow
//==============================================================================

void applyAtmosphere( inout vec3 col, in CRay r, in CHit hit )
{
#ifdef ENABLE_FOG
    float fogAmount = exp(hit.t * -fogDensity);
    vec3 fogCol = getSkyGradient(r.d);
    #ifdef ENABLE_DIRECTIONAL_LIGHT_FLARE
        CDirLight dirLight = getDirLight();
        float dirDot = clamp(dot(-dirLight.dir, r.d), 0.0, 1.0);
        fogCol += dirLight.col * pow(dirDot, 10.0);
    #endif
    col = mix(fogCol, col, fogAmount);
#endif

#ifdef ENABLE_POINT_LIGHT_FLARE
    CPointLight pointLight = getPointLight();
    vec3 toLight = pointLight.pos - r.o;
    float pointDot = dot(toLight, r.d);
    pointDot = clamp(pointDot, 0.0, hit.t);
    vec3 PointClosest = r.o + r.d * pointDot;
    float dist = length(PointClosest - pointLight.pos);
    col += pointLight.col * 0.01/ (dist * dist);
#endif
}

float Schlick( in vec3 h, in vec3 v, in float ior, in float smoothFactor )
{
	float F0 = pow( (airIOR - ior)/(airIOR + ior), 2.0 ); // airIOR = 1.0003
	float dotHV = dot(h, v); // H o V
    return F0 + (1.0 - F0) * pow(saturate(1.0 - dotHV), 5.0) * smoothFactor;// fresnel = F0 + (1.0 - F0) * (1 - cosTheta)^5 * smoothFactor
}

vec3 getFresnelColor( in vec3 diffuse, in vec3 specular, in vec3 n, in vec3 rd, in CMaterial m )
{
	vec3 ld = reflect(rd, n); // R = reflect(-V, N)
	vec3 h = normalize(ld - rd); // H = L + V
	//float F = Schlick(h, -rd, m.ior, m.smoothness * 0.9 + 0.1);
	float F = Schlick(h, -rd, m.ior, m.smoothness);
	return mix(diffuse, specular, F); // return (Ld*Md) * (1-F) + Ls * (F)
}

float getSpecularIntensity( in vec3 rd, in vec3 ld, in vec3 n, in float smoothness )
// ld = surface to light position
{          
	vec3 h = normalize(ld - rd);  // H = L + V
	float dotNH = max(0.0, dot(n, h));      // NoH = N o H
	float specPower = exp2(4.0 + 6.0 * smoothness);   // smoothness ==> specPower (= shininess)
	float specIntensity = (specPower + 2.0) * 0.125; // specPower ==> specIntensity (= specIntensity)
	return pow(dotNH, specPower) * specIntensity;   // specWt = specIntensity * NoH^(shininess)
}

CShading getPointLighting( in CPointLight light, in vec3 p, in vec3 n, in vec3 rd, in CMaterial m )
{
	vec3 ld = normalize(light.pos - p);
	float tmin = 0.01 / abs(dot(ld, n));
	float tmax = length(light.pos - p);
	#ifdef ENABLE_SHADOW
		float shadow = sceneShadow( p, ld, tmin, tmax, 32.0 );
	#else
		float shadow = 1.0;
	#endif

	CShading shading;
	float atten = 1.0 / (tmax * tmax); // atten = 1/(d*d)

	vec3 inLight = light.col * (shadow*atten) * max(0.0, dot(ld, n)); // Ld = shadow * atten * Lc * NoL
	shading.diffuse = inLight;
	shading.specular = getSpecularIntensity( rd, ld, n, m.smoothness ) * inLight; // Ls = Ld * specWt
	return shading;
}

CShading getDirectionalLighting( in CDirLight light, in vec3 p, in vec3 n, in vec3 rd, in CMaterial m )
{
	vec3 ld = light.dir;
	float tmin = 0.01 / abs(dot(ld, n));
	float tmax = FAR*0.5;
	#ifdef ENABLE_SHADOW
		float shadow = sceneShadow( p, ld, tmin, tmax, 32.0 );
	#else
		float shadow = 1.0;
	#endif

	CShading shading;
	vec3 inLight = light.col * (shadow) * max(0.0, dot(ld, n)); // Ld = shadow * Lc * NoL
	shading.diffuse = inLight;
	shading.specular = getSpecularIntensity( rd, ld, n, m.smoothness ) * inLight; // Ls = Ld * specWt
	return shading;
}

float getSubSurfScattering( in vec3 p, in vec3 rd )
// subsurface scattering
{
	float SSS_K = 1.5;
	float SSS_DELTA = 0.3;
	const int SSS_N = 5;
	float sum = 0.0;
	float weight = -0.5;
	float delta = SSS_DELTA;
	for(int i = 0; i < SSS_N; i++)
	{
		sum += weight * min( 0.0, sceneMap(p + delta*rd).x );
		delta += delta;
		weight *= 0.5;
	}
	return clamp( float(SSS_K)*sum, 0.0, 1.0 );
}

vec3 getPhongShading( in CRay r, in CHit hit, in CMaterial m )
// r = 1st ray (incident)
{
	#ifdef ENABLE_BUMP_MAPPING
		if( hit.m == kMaterialTexture0 )
			hit.n = getBumpNormal(textureMaps[0], hit.p, hit.n, 0.075);
		else if( hit.m == kMaterialTexture1 )
			hit.n = getBumpNormal(textureMaps[1], hit.p, hit.n, 0.075);
		else if( hit.m == kMaterialTexture2 )
			hit.n = getBumpNormal(textureMaps[2], hit.p, hit.n, 0.075);
		else if( hit.m == kMaterialTexture3 )
			hit.n = getBumpNormal(textureMaps[3], hit.p, hit.n, 0.075);
	#endif

	vec3 sceneCol;
	CShading shading;
	shading.diffuse = vec3(0.0);  // totalLd = 0.0
	shading.specular = vec3(0.0); // totalLs = 0.0

	#ifdef ENABLE_AO
		float ao = sceneAO(hit.p, hit.n);
	#else
		float ao = 1.0;
	#endif

	// reflection color added to specular term
	vec3 ambientLight = getAmbientLight(hit.n);
	shading.diffuse += ambientLight * ao;
	shading.specular += hit.refCol;

	#ifdef ENABLE_SUNDIR_LIGHT
		CDirLight sunLight = getSunLight();
		CShading sunShading = getDirectionalLighting(sunLight, hit.p, hit.n, r.d, m);
		shading.diffuse += sunShading.diffuse * ao;
		shading.specular += sunShading.specular;
	#endif

	#ifdef ENABLE_DIRECTIONAL_LIGHT
		CDirLight dirLight = getDirLight();
		CShading dirShading = getDirectionalLighting(dirLight, hit.p, hit.n, r.d, m);
		shading.diffuse += dirShading.diffuse * ao;
		shading.specular += dirShading.specular;
	#endif

	#ifdef ENABLE_POINT_LIGHT
		CPointLight pointLight = getPointLight();
		CShading pointShading = getPointLighting(pointLight, hit.p, hit.n, r.d, m);
		shading.diffuse += pointShading.diffuse * ao;
		shading.specular += pointShading.specular;
	#endif

	vec3 diffuseReflected = shading.diffuse * m.albedo;
	diffuseReflected = mix(diffuseReflected, hit.trmCol, m.transparency);
	#ifdef ENABLE_SPECULAR
		sceneCol = getFresnelColor(diffuseReflected, shading.specular, hit.n, r.d, m);
	#else
		sceneCol = diffuseReflected;
	#endif

	#ifdef ENABLE_SUBSURFACE
		float sss = getSubSurfScattering( hit.p, r.d );
		sceneCol *= (1.0 - sss);
	#endif

	return sceneCol;
}

vec3 getSceneColor2( in CRay r );

vec3 getReflection( in CRay r, in CHit hit )
// r = 1st ray (incident)
{
#if defined (ENABLE_SPECULAR) && defined (ENABLE_REFLECTIONS)
    float separation = 0.01;
    CRay re; // 2nd ray (reflected)
    re.d = reflect(r.d, hit.n);
    re.o = hit.p;
	re.tmin = separation / abs(dot(re.d, hit.n));
    re.tmax = FAR*0.5;
    return getSceneColor2(re);
#else
    return getSkyGradient(reflect(r.d, hit.n));                              
#endif
}

vec3 getTransmission( in CRay r, in CHit hit, in CMaterial m )
// r = 1st ray (incident)
{
	#ifdef ENABLE_TRANSPARENCY
	{
		// Trace until outside transparent object
		float separation = 0.01;
		CRay ra;
		ra.d = refract(r.d, hit.n, airIOR/m.ior);
		ra.o = hit.p;
		ra.tmin = separation / abs(dot(ra.d, hit.n));
		ra.tmax = FAR*0.5;

		#ifdef ENABLE_DOUBLE_TRANSPARENCY
			CHit hit2;
			rayMarching(ra, hit2);
			vec3 n = sceneNormal(hit2.p);

			CRay ra2;
			ra2.d = refract(ra.d, n, m.ior/airIOR);
			ra2.o = hit2.p;
			ra2.tmin = separation / abs(dot(ra2.d, n));
			ra2.tmax = FAR*0.5;
			
			float extinctionDist = hit2.t;
			vec3 sceneCol = getSceneColor2(ra2);
		#else
			float extinctionDist = 0.5;
			vec3 sceneCol = getSceneColor2(ra);
		#endif

		vec3 extinctionCol = m.albedo;
		// extinction should really be exp(-) but this is a nice hack to get RGB
		extinctionCol = (1.0 / (1.0 + (extinctionCol * extinctionDist)));
		return sceneCol * extinctionCol;
	}
	#else
        return getSkyGradient(reflect(r.d, hit.n));
    #endif
}

vec3 getSceneColor2( in CRay r )
// no reflections, no transparency
// r = secondary rays
{
	CHit hit;
	rayMarching(r, hit);

	vec3 sceneCol;
	if( hit.m == kMaterialSky )
	{
		sceneCol = getSkyGradient(r.d);
	}
	else
	{
		hit.n = sceneNormal(hit.p);
		CMaterial m = getMaterial(hit);

		// use sky gradient instead of reflection
		hit.refCol = getSkyGradient(reflect(r.d, hit.n));
		m.transparency = 0.0;

		sceneCol = getPhongShading(r, hit, m);
	}

	applyAtmosphere(sceneCol, r, hit);
	return sceneCol;
}

vec3 getSceneColourTestVersion( in CRay r )
{
	CHit hit;
	rayMarching(r, hit);

	vec3 sceneCol;
	if( hit.m == kMaterialSky )
	{
		sceneCol = getSkyGradient(r.d);
	}
	else
    {
		hit.n = sceneNormal(hit.p);

		#ifdef ENABLE_AO
			float ao = sceneAO(hit.p, hit.n);
		#else
			float ao = 1.0;
		#endif

		#ifdef ENABLE_SHADOW
			float shadow = sceneShadow(hit.p, SUN_LIGHT);
		#else
			float shadow = 1.0;
		#endif

		//vec3 ambientLight = getAmbientLight(hit.n);

		CMaterial m = getMaterial(hit);
		sceneCol = getPhongShading( m.albedo, 32.0, SUN_LIGHT, SUN_COLOR, hit.p, hit.n, r.d );
		sceneCol *= (shadow * ao);
    }
	return sceneCol;
}

vec3 getSceneColor( in CRay r )
{
#ifdef ENABLE_TEST_RENDER
	return getSceneColourTestVersion( r );
#endif

	// step 1: geometric intersection
	CHit hit;
	rayMarching(r, hit);

	vec3 sceneCol;
	if( hit.m == kMaterialSky )
	{
		// no hit
		sceneCol = getSkyGradient(r.d);
		gl_FragDepth = 0.99;
	}
	else
	{
		// step 2: normal
		hit.n = sceneNormal(hit.p);

		// step 3: material
		CMaterial m = getMaterial(hit);

		// step 4.1: 2nd ray (reflected)
		hit.refCol = getReflection(r, hit);

		// step 4.2: 2nd ray (transmitted)
		hit.trmCol = (m.transparency > 0.0)? getTransmission(r, hit, m) : vec3(0.0);

		// step 4.3: 1st ray (incident)
		sceneCol = getPhongShading(r, hit, m);

		gl_FragDepth = getFragDepth( hit.p );
	}

	// step 5: fog...
	applyAtmosphere(sceneCol, r, hit);
	return sceneCol;
}

void applyPostEffects( inout vec3 col )
{
	// tonemap
	col = (1.0 - exp(-col * 1.5)) * 1.0024; // fExposure = 1.5

	// gamma correction
	col = pow( col, vec3(0.4545) );
}

//==============================================================================
// Main
//==============================================================================

void main()
{
	@import ./ray;

	CRay r;
	r.o = ro;
	r.d = rd;

#if defined (DEMO_REPEAT_SPHERES)
	r.tmin = 0.0;
	r.tmax = FAR;
#elif defined (DEMO_PRIMITIVES)
	r.tmin = 0.0;
	r.tmax = 20.0;
#endif

	vec3 col = getSceneColor( r );
	applyPostEffects( col );
	gl_FragColor = vec4( col, 1.0 );
}