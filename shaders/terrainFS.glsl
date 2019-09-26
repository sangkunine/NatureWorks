@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;
#define USE_DUSTWIND_1
#define USE_TERRAIN_1
@import ./scene;

// textureMaps[0] = 'images/raymarch/grayNoise256.png'
#define iChannel0   textureMaps[0]

// on the derivatives based noise: http://iquilezles.org/www/articles/morenoise/morenoise.htm
// on the soft shadow technique: http://iquilezles.org/www/articles/rmshadows/rmshadows.htm
// on the fog calculations: http://iquilezles.org/www/articles/fog/fog.htm
// on the lighting: http://iquilezles.org/www/articles/outdoorslighting/outdoorslighting.htm
// on the raymarching: http://iquilezles.org/www/articles/terrainmarching/terrainmarching.htm

// Choose...
#define ENABLE_AUTO_VIEW
#define ENABLE_DUST_WIND

const float FAR = 300.0;//200.0
const vec3 SUN_LIGHT = normalize(vec3(-0.8, 0.4, -0.3));
const vec3 SUN_COLOR = vec3(1.0, 0.9, 0.85);
const float SKY_HEIGHT = 300.0;

float skyMap( in vec3 p )
{
    return SKY_HEIGHT - p.y;
}

vec2 sceneMap( in vec3 p )
{
    vec2 res = vec2( skyMap(p), MATERIAL_SKY );
    res = dUnion( res, vec2( terrainMap( iChannel0, p ), MATERIAL_TERRAIN ) );
    return res;
}

void rayMinMax( in vec3 ro, in vec3 rd, out float tmin, out float tmax )
{
    tmin = 1.0;
    tmax = FAR;
    float max_height = SKY_HEIGHT;
    float t = (max_height - ro.y) / rd.y;
    if( t > 0.0 )
    {
        if( ro.y > max_height ) tmin = max( tmin, t );
        else                    tmax = min( tmax, t );
    }
    else
    {
        if( ro.y > max_height ) tmin = tmax = 1.0;
    }
}

vec3 render( in vec3 ro, in vec3 rd )
{
    float tmin, tmax;
	rayMinMax( ro, rd, tmin, tmax );

    vec2 tm = bisectMarching( ro, rd, tmin, tmax );
    if( tm.x > tmax ) tm.y = MATERIAL_SKY;

    vec3 col;
    if( tm.y == MATERIAL_SKY )
    {
        gl_FragDepth = 0.99;

        col = skyColor( SUN_LIGHT, rd, 0.0 );//1.0
        applyClouds( col, iChannel0, ro, rd );
	}
    else if( tm.y == MATERIAL_TERRAIN )
    {
        vec3 p = ro + rd * tm.x;
        gl_FragDepth = getFragDepth( p );

        col = terrainColor( iChannel0, SUN_COLOR, SUN_LIGHT, p, tm.x, tmax );
        // height-based fog density
        applyFog( col, SUN_COLOR, SUN_LIGHT, ro, rd, tm.x*1.4 );
        // constant fog density
        applyFog( col, SUN_COLOR, SUN_LIGHT, rd, 0.003, tm.x );//0.005
    }

    // sun scatter
    col += sunScatter( SUN_LIGHT, rd );

#ifdef ENABLE_DUST_WIND
    float dustAmount = 0.25;//0.25;
    float dustHeight = 75.0*TERRAIN_SCALE;//50.0
    float windTurbulency = 0.5;
    applyDustWind( col, ro, rd, tm.x, dustAmount, dustHeight, windTurbulency );
#endif

    col = FilmicToneMapping( col );
    //col = col * 1.05 - 0.05;
    col = LinearToGamma( vec4(col, 1.0), 0.8 ).rgb;//0.8

    return col;
}

#ifdef ENABLE_AUTO_VIEW
void cameraAutoView( in sampler2D tex, out vec3 ro, out vec3 rd )
{
    float curTime = 5.5*time;
    ro = vec3( 0.0, 0.0, -95.0-curTime );
    vec3 ta = vec3( 0.0, 0.0, -110.0-curTime );
    ta = mix( ro + vec3(0.0, 1.0, 0.0), ta, smoothstep(1.0, 25.0, curTime) );
    ro.y = terrainL( tex, ro.xz ) + 50.0*TERRAIN_SCALE;//75.0
    ta.y = ro.y - 40.0*TERRAIN_SCALE;//-75.0

    float fl = 1.0;
    vec2 xy = (-resolution.xy + 2.0*gl_FragCoord.xy)/resolution.y;
    mat3 cam = cameraMatrix( ro, ta ); // cam[0] = cu, cam[1] = cv, cam[2] = cw
    rd = normalize( xy.x*cam[0] + xy.y*cam[1] + fl*cam[2] );
}
#endif

void main()
{
#ifdef ENABLE_AUTO_VIEW
    vec3 ro, rd;
    cameraAutoView( iChannel0, ro, rd );
#else
    @import ./ray;
#endif

    vec3 col = render( ro, rd );
    gl_FragColor = vec4( col, 1.0 );
    gl_FragColor.rgb = Vignetting( gl_FragColor.rgb, 0.5 );
}