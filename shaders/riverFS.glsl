@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;

#define NEW_RIVER_DEPTH
const float RIVER_WAVE_LENGTH = 20.0;
const float RIVER_WAVE_HEIGHT = 1.2;
float riverCurve( float x ) {
    // definition of meandering river
    float L = RIVER_WAVE_LENGTH;
    float H = RIVER_WAVE_HEIGHT;
    //return H * sin( PI2/L * x );
    float W = PI2/L;
    return H*(sin(W*x) + sin(W*x*2.0));
}
float riverCurved( float x ) {
    // derivative of riverCurve() above
    float L = RIVER_WAVE_LENGTH;
    float H = RIVER_WAVE_HEIGHT;
    float W = PI2/L;
    //return H * W * cos(W * x);
    return H*W*(cos(W*x) + 2.0*cos(W*x*2.0));
}
float riverBaseDepth( vec3 p ) {
    float depth = 0.3 + (0.5 + 0.5*sin(p.x*0.001 + 3.0)) * 0.4;
    float width = 2.0 + cos( p.x * 0.1 ) * 1.0;
    float amount = smoothstep( width, width * 0.5, abs(p.z - riverCurve(p.x)) );
    return amount * depth;
}
float riverDepth( vec3 p ) {
    float depth = 0.0;
    depth += riverBaseDepth(p);
    depth += 0.5*riverBaseDepth(p*2.0);
    //depth += 0.25*riverBaseDepth(p*4.0);
    return depth;
}
#define USE_WATER_2
@import ./scene;

// textureMaps[0] = 'images/raymarch/grayNoise256.png'
// textureMaps[1] = 'images/raymarch/lichen.jpg'
#define iChannel0   textureMaps[0]
#define iChannel1   textureMaps[1]

//==============================================================================

// Choose...
#define ENABLE_AUTO_VIEW

const float FAR = 60.0;//20.0;
const vec3 SUN_LIGHT = normalize(vec3(-1.0, 0.7, 0.25));
const vec3 SUN_COLOR = vec3(1.0, 0.85, 0.5) * 5.0;//5.0;
const vec3 SKY_COLOR = vec3(0.1, 0.5, 1.0);

float terrainHeight( vec3 p )
// p.y not used...
{
    // terrain map
    vec2 q = p.xz * vec2(0.5, 1.0);
    float f = 0.0;
    float tot = 0.0;
    float a = 1.0;
    for(int i=0; i<2; i++)//3
    {
        f += noise( q ) * a;
        q *= 2.0;
        tot += a;
        a *= 0.5;
    }
    float fbm = f / tot;
    float height = fbm * fbm;

    // add river...
    height -= riverDepth(p);

    return height;
}

vec2 sceneMap( in vec3 p ) // mandatory function
{
    return vec2( p.y - terrainHeight( p ), MATERIAL_TEXTURE1 );
}

vec3 getTerrainColor( in vec3 lc, in vec3 ld, in vec3 sky, in vec3 ro, in vec3 rd, in float t )
{
    // position & normal
    vec3 p = ro + rd * t;
    vec3 n = sceneNormal( p, EPS );
    gl_FragDepth = getFragDepth( p );

    // material
    vec3 col = texCube( iChannel1, p, n ).xyz;
    vec3 albedo = col*col;
    float wetness = 1.0 - saturate( (p.y + 0.025) * 5.0 );
    float gloss = mix( albedo.r, 1.0, wetness );
    albedo = mix( albedo, albedo * 0.8, wetness );
    vec3 specF0 = vec3(0.001);

    // sun light ==> diffuseCol + specularCol
    vec3 diffuseCol;
    vec3 specularCol;
    float shadow = sceneShadow( p, SUN_LIGHT, FAR*0.01, FAR*0.1, 16.0 );
    getCookShading( albedo, gloss, SUN_COLOR, SUN_LIGHT, n, rd, diffuseCol, specularCol );
    diffuseCol *= shadow;

    // reflectCol ==> specularCol
    vec3 reflectRd = reflect( rd, n );
    vec3 reflectCol = reflectRayColor( lc, ld, rd, p, n, FAR ).xyz;
    reflectCol = mix( waterEnvColor(sky, reflectRd, gloss), reflectCol, pow(gloss, 40.0) );
    specularCol += reflectCol;

    // fresnel
    vec3 fresnel = fresnelGloss( n, -rd, specF0, gloss );
    float specScale = 1.0;

    vec3 result = mix( diffuseCol, specularCol, fresnel*specScale );
    return result;
}

vec3 sceneColor( in vec3 ro,  in vec3 rd )
{
    vec2 hit = bisectMarching( ro, rd, 0.01, FAR );

    float fogDistance;
    vec3 result;

    // sky
    if( hit.x > FAR )
    {
        result = skyColor( SUN_LIGHT, rd, 0.0 );
        fogDistance = FAR;
        gl_FragDepth = 0.99;
    }

    // water
    else
    {
        vec4 waterCol = waterColor( SUN_COLOR, SUN_LIGHT, SKY_COLOR, ro, rd, hit.x, FAR );
        result = waterCol.xyz;
        fogDistance = waterCol.w;
    }

    // terrain
    if( result == vec3(0.0) )
    {
        result = getTerrainColor( SUN_COLOR, SUN_LIGHT, SKY_COLOR, ro, rd, hit.x );
        fogDistance = hit.x;
    }

    // fog...
    float fogDensity = 0.025;//0.015
    vec3 fogCol = skyColor( SUN_LIGHT, rd, 0.0 );
    applyFog( result, fogCol, fogDensity, fogDistance );

    return result;
}

#ifdef ENABLE_AUTO_VIEW
mat3 cameraAutoView( out vec3 ro, out vec3 rd )
{
    vec3 ta = vec3(0.0, -0.5, 0.0);
    ta.x -= time * 0.75;//0.5
    ro = ta + vec3(1.5);//2.0
    
    float fHeading = time * 0.15;//0.1
    float dist = 1.5 - cos(time * 0.1 + 2.0) * 0.8;

    ro.x += sin( fHeading ) * dist;
    ro.z += cos( fHeading ) * dist;
    ro.y += 1.0 + dist * dist * 0.01;
    ta.z += riverCurve( ta.x );
    ro.z += riverCurve( ro.x );
#if 1
    // camera moves above the water
    ro.y = max( ro.y, terrainHeight( ro ) + 0.2 );
#else
    // camera moves up and down the water
    ro.y = min( ro.y, terrainHeight( ro ) + 0.5 );
#endif

    // camera ray direction
    float fl = 1.5;//2.0
    vec2 xy = (-resolution.xy + 2.0*gl_FragCoord.xy)/resolution.y;
    mat3 cam = cameraMatrix( ro, ta ); // cam[0] = cu, cam[1] = cv, cam[2] = cw
    rd = normalize( xy.x*cam[0] + xy.y*cam[1] + fl*cam[2] );
    return cam;
}
#endif

void main()
{
#ifdef ENABLE_AUTO_VIEW
	vec3 ro, rd;
    cameraAutoView( ro, rd );

#else
    @import ./ray;

#endif

    vec3 result = sceneColor(ro, rd);
	result = Vignetting( result, 0.5 );
    result = FilmicToneMapping( result );
    result = LinearToGamma( vec4(result, 1.0), 0.9 ).xyz; //2.2

	gl_FragColor = vec4(result, 1.0);
}