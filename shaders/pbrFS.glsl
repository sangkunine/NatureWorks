//==============================================================================
//
// Microsurface (roughness or glossiness)
// -metalness workflow: setting a metalness value or creating a metalness map(1 for metal, 0 for non-metal)
// -specular workflow: setting a specular value or creating a specular map(color and intensity of the reflected light)
//
// Diffuse(albedo or base color): 
// -constant value (= Gold(1.000, 0.766, 0.336), Silver(0.972, 0.960, 0.915), Copper(0.955, 0.637, 0.538)...)
// -diffuse map texture
//
// Metalness:
// -single value (metal (1) or non-metal(0))
// -metalness map (= define specific areas of your material as metal or non-metal)
//
// Glossiness:
// -blurry or sharp the reflections on the material (or how broad or narrow the specular highlights are)
// -single value between (0-100) or a glossiness map
// (*The roughness is the inverse of the glossiness)
//
// Additional properties: Ambient Occulusion, Emissive, Opacity, Normal and Height maps
//
//==============================================================================

@import ./include;
@import ./utils;
@import ./distances;

//==============================================================================

//#undef STANDARD
#define USE_RAY_SHADOW
#define USE_RAY_AO

@import ./pbr;

//==============================================================================

// numMaterials: must be <= 4
#define numMaterials    4

// Choose...
#define ENABLE_BUMP_NORMAL
#define ENABLE_FOG
#define ENABLE_REFLECT_PASS

// Choose one of these demos...
#define ENABLE_DEMO_ROUGH_METAL
//#define ENABLE_DEMO_PBR_MATERIALS
//#define ENABLE_DEMO_PRIMITIVES

// floor material
#define MTL_FLOOR       100.0
#define FLOOR_COLOR     checkerColor(p, 5.0, 0.6)

#ifdef ENABLE_FOG
const float fogDensity = 0.001;
const vec3 fogColor = vec3(0.56, 0.64, 0.8);
#endif

// const vec3 cGold = vec3(1.000, 0.766, 0.336);
// const vec3 cSilver = vec3(0.972, 0.960, 0.915);
// const vec3 cCopper = vec3(0.955, 0.637, 0.538);
// const vec3 cIron = vec3(0.56, 0.57, 0.58);
// const vec3 cAluminium = vec3(0.913, 0.921, 0.925);
// const vec3 Chromium = vec3(0.550, 0.556, 0.554);
// const vec3 Nickel = vec3(0.660, 0.609, 0.526);
// const vec3 Titanium = vec3(0.542, 0.497, 0.449);
// const vec3 Cobalt = vec3(0.662, 0.655, 0.634);
// const vec3 Platinum = vec3(0.672, 0.637, 0.585);


//==============================================================================
// Scene Modeling
//==============================================================================

vec2 simpleScene(vec3 p)
{
    // // one sphere
    // vec2 res = vec2(fPlane(p), MTL_FLOOR);
    // float r = 0.5;
    // vec3 c = vec3(0.0, r, 0.0);
    // vec2 sph = vec2(fSphere(p - c, r), 0.0);
    // res = dUnion(res, sph);
    // return res;

    // smooth union
    vec2 res = vec2(fPlane(p), MTL_FLOOR);
    float k = 0.1;
	vec3 tp = sTranslate( p, vec3(0.0, 2.0, 0.0) );
	vec2 dcy1 = vec2( fCylinder(tp, vec2(1.0, 2.0)), 0.0 );
	vec2 dcy2 = vec2( fCylinder(tp.yzx, vec2(0.5, 2.0)), 1.0 );
	//res = dUnion( res, dSmoothUnion( dcy1, dcy2, k ) );
    res = dUnion( res, dSmoothSubtract( dcy1, dcy2, k ) );
    return res;
}

#ifdef ENABLE_DEMO_ROUGH_METAL
vec2 demoRoughMetalScene(vec3 p)
{
    vec2 res = vec2( fPlane(p), MTL_FLOOR ); // floor

    const int nx = 7, ny = 7;
    float r = 0.25;
    float m = 0.0;

    float dx = 3.0*r, dy = 3.0*r;
    float x, x0 = -(float(nx)/2.0-0.5)*dx;
    float y, y0 = r;

    for(int j = 0; j < ny; j++) {
        y = y0 + float(j)*dy;
        for(int i = 0; i < nx; i++) {
            x = x0 + float(i)*dx;
            res = dUnion( res, vec2( fSphere(p-vec3(x, y, 0.0), r), m ) );
            m = m + 1.0;
        }
    }
    return res;
}
#endif

#ifdef ENABLE_DEMO_PBR_MATERIALS
vec2 demoPBRMaterialsScene(vec3 p)
{
    vec2 res = vec2(fPlane(p), MTL_FLOOR);
    float r = 0.5;
    //float dr = 0.0;
    float dr = 0.2;
    float dx = (r + dr)*2.0;
    float x0 = (float(numMaterials)-1.0) / -2.0;
    for(int i = 0; i < numMaterials; i++) {
        vec3 c = vec3( x0 + float(i)*dx, r, 0.0 );
        vec2 sph = vec2(fSphere(p - c, r), float(i));
        res = dUnion(res, sph);
    }
    return res;
}
#endif

#ifdef ENABLE_DEMO_PRIMITIVES
vec2 demoPrimitivesScene(vec3 p)
{
    // 21 objects + 1 floor
    float m = 0.0;
    vec2 res = dUnion( vec2( fPlane(     p), MTL_FLOOR ),
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

vec2 sceneMap(vec3 p)
{
#if defined(ENABLE_DEMO_ROUGH_METAL)
    return demoRoughMetalScene(p);

#elif defined(ENABLE_DEMO_PBR_MATERIALS)
    return demoPBRMaterialsScene(p);

#elif defined(ENABLE_DEMO_PRIMITIVES)
    return demoPrimitivesScene(p);

#else
    return simpleScene(p);

#endif
}

//==============================================================================
// Material
//==============================================================================

void getPBRMaterial( in float m, in vec3 p, in vec3 n, out PhysicalMaterial material )
{
    //------------------------------------------------------
    vec3 baseColor = vec3(1.0);         // default baseColor
    float metalness = 0.1;//1.0;           // default metalness
    float roughness = 0.1;//0.04;  // default roughness
    //------------------------------------------------------

    if( m == MTL_FLOOR )
    {
        baseColor = FLOOR_COLOR;
        metalness = 0.1;
        roughness = 0.1;
    }
    else
    {
//------------------------------------------------------
#if defined(ENABLE_DEMO_ROUGH_METAL)
        baseColor = vec3(1.0, 1.0, 1.0);
        float nx = 7.0, ny = 7.0;
        metalness = floor(m/nx) / (ny-1.0); // metalness = [0,1] along yaxis
        roughness = mod(m,nx) / (nx-1.0);   // roughness = [0,1] along xaxis

//------------------------------------------------------
#elif defined(ENABLE_DEMO_PRIMITIVES)
        baseColor = hsl2rgb( vec3( m/21.0, 1.0, 0.5) );
        metalness = 0.7;//m/21.0;
        roughness = n.y;

//------------------------------------------------------
#else
        for(int i = 0; i<numMaterials; i++)
        {
            if( int(m) == i ) {
                vec3 texCol = pow( texCube(albedoMaps[i], p, n).rgb, vec3(2.2) );
                baseColor = mix( albedos[i], texCol, float(useAlbedoMaps[i]) );
                metalness = mix( metalnesss[i],  texCube(metalnessMaps[i], p, n).r, float(useMetalnessMaps[i]) );
                roughness = mix( roughnesss[i], texCube(roughnessMaps[i], p, n).r, float(useRoughnessMaps[i]) );
                break;
            }
        }
#endif
//------------------------------------------------------
    }

    // material.metalness
    material.metalness = clamp( metalness, 0.0, 1.0 );

    // material.specularRoughness <== roughness (0.04 ~ 1.0)
    material.specularRoughness = clamp( roughness, 0.04, 1.0 );

    // material.diffuseColor <== baseColor + metalness
    material.diffuseColor = baseColor.rgb * ( 1.0 - material.metalness );

    // material.specularColor <== baseColor + metalness
    #ifdef STANDARD
        material.specularColor = mix( vec3( DEFAULT_SPECULAR_COEFFICIENT ), baseColor.rgb, material.metalness );
    #else
        float reflectivity = 0.5; // maps to F0 = 0.04
        material.specularColor = mix( vec3( MAXIMUM_SPECULAR_COEFFICIENT * pow2( reflectivity ) ), baseColor.rgb, material.metalness );
        material.clearCoat = saturate( clearCoat ); // Burley clearcoat model
        material.clearCoatRoughness = clamp( clearCoatRoughness, 0.04, 1.0 );
    #endif
}

//==============================================================================
// Post-Effects
//==============================================================================

#ifdef ENABLE_BUMP_NORMAL
vec3 getBumpNormal( in float m, in vec3 p, in vec3 n, float bumpFactor )
// bumpFactor = 0.0075, 0.075
{
    vec3 normal = n;
    for(int i = 0; i<numMaterials; i++)
    {
        if( int(m) == i )
        {
            if( useAlbedoMaps[i] == 1 )
                normal = getBumpNormal( albedoMaps[i], p, n, bumpFactor );
            break;
        }
    }
    return normal;
}
#endif

#ifdef ENABLE_FOG
vec3 applyFog( in vec3 color, in vec3 fogColor, in float fogDensity, in float fogDepth )
{
    float fogFactor = exp2( -fogDensity * fogDensity * fogDepth * fogDepth * LOG2 );//LOG2 = 1.442695
    fogFactor = 1.0 - clamp( fogFactor, 0.0, 1.0 );
	return mix( color, fogColor, fogFactor );
}
#endif

void doPostEffects( inout vec4 color, in float fogDepth )
{
    color.rgb = FilmicToneMapping( color.rgb );

#ifdef ENABLE_FOG
    color.rgb = applyFog( color.rgb, fogColor, fogDensity, fogDepth );
#endif

    color = LinearToGamma( color, 2.2 ); // gammaFactor = 2.2

    color.rgb *= color.a;

    color.rgb = Dithering( color.rgb );
}

//==============================================================================
// Rendering
//==============================================================================

vec4 render(in vec3 ro, in vec3 rd, out float t)
{
    vec4 color = vec4(0.0);

    // hit.xy = (distance, material)
    vec2 hit = rayMarching(ro, rd);
    t = hit.x;
    float m = hit.y;

    gl_FragDepth = 1.0;

    // hit something in the scene
    if (m >= 0.0)
    {
        vec3 p = ro + t * rd;
        vec3 n = sceneNormal(p);
        gl_FragDepth = getFragDepth( p );

        PhysicalMaterial material;
        getPBRMaterial( m, p, n, material );

        #ifdef ENABLE_BUMP_NORMAL
            n = getBumpNormal( m, p, n, material.specularRoughness * 0.01 );
        #endif

        color = getPBRColor( p, n, -rd, material );

        #ifdef ENABLE_REFLECT_PASS
        if( material.metalness > 0.0 )
        {
            vec3 r = reflect(rd, n);
            vec2 hit2 = rayMarching(p, r);
            m = hit2.y;
            if( m >= 0.0 )
            {
                p += r * hit2.x;
                n = sceneNormal(p);

                PhysicalMaterial material2;
                getPBRMaterial( m, p, n, material2 );

                float gloss = material.metalness;
                gloss *= (1.0 - material.specularRoughness);
                gloss *= 0.35 / (1.0 + hit2.x * hit2.x);
                color += getPBRColor( p, n, -r, material2 ) * gloss;
            }
        }
        #endif
    }

    return color;
}

//==============================================================================
// Main
//==============================================================================

void main()
{
    @import ./ray;

    float t;
    vec4 color = render( ro, rd, t );

    if( gl_FragDepth == 1.0 ) discard;

    doPostEffects( color, t );

    gl_FragColor = color;
}