@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;
#define USE_CLOUDS_2
#define USE_TERRAIN_2
#define USE_TREES_1
@import ./scene;

// textureMaps[0] = 'images/raymarch/grayNoise256.png'
#define iChannel0   textureMaps[0]

// Choose...
#define ENABLE_AUTO_VIEW

const float FAR = 500.0;//300.0;
const vec3 SUN_LIGHT = vec3(-0.624695,0.468521,-0.624695);
const vec3 SUN_COLOR = vec3(1.0, 0.9, 0.85);
const float SKY_HEIGHT = 500.0;

float skyMap( in vec3 p )
{
    return SKY_HEIGHT - p.y;
}

vec2 sceneMap( in vec3 p )
{
    vec2 res = vec2( skyMap(p), MATERIAL_SKY );
    res = dUnion( res, vec2( terrainMap( iChannel0, p ), MATERIAL_TERRAIN ) );
    res = dUnion( res, vec2( treesMap( iChannel0, p ), MATERIAL_TREES ) );
    //res = dUnion( res, vec2( cloudsMap( p ).x, MATERIAL_CLOUDS ) );
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

    // raymarch ==> intersection
    vec2 tm = bisectMarching( ro, rd, tmin, tmax );
    if( tm.x > tmax ) tm.y = MATERIAL_SKY;

    vec3 col = vec3(0.0);
    if( tm.y == MATERIAL_SKY )
    {
        gl_FragDepth = 0.99;

        // sky & clouds
        col = skyCloudsColor( iChannel0, SUN_COLOR, SUN_LIGHT, ro, rd );

        // clouds(volumetric)
        vec4 cloudsCol = cloudsColor( SUN_LIGHT, ro, rd, 0.0, FAR );//1000.0
        #if 0
            col = col*(1.0-cloudsCol.w) + cloudsCol.xyz;
        #else
            col = mix( col, cloudsCol.xyz, cloudsCol.w );
        #endif
	}

    else if( tm.y == MATERIAL_TERRAIN )
    {
        vec3 p = ro + rd * tm.x;
        gl_FragDepth = getFragDepth( p );
        col = terrainColor( iChannel0, SUN_LIGHT, rd, p, tm.x );
    }
    else if( tm.y == MATERIAL_TREES )
    {
        vec3 p = ro + rd * tm.x;
        gl_FragDepth = getFragDepth( p );
        col = treesColor( iChannel0, SUN_LIGHT, p, rd, tm.x );
    }

    // fog
    applyFog( col, FOG_COLOR, 0.003, tm.x*0.5 );

    // sun glare
    col += sunScatter( SUN_LIGHT, rd );

    // color grading
#if 0
    col = col*0.15 + 0.85*col*col*(3.0-2.0*col); // contrast
    col = pow( col, vec3(1.0,0.92,1.0) );  // soft green
    col *= vec3(1.02,0.99,0.99);           // tint red
    col.z = (col.z+0.1)/1.1;               // bias blue
    col = mix( col, col.yyy, 0.15 );       // desaturate
    col = saturate( col );
#endif

    col = FilmicToneMapping( col );
    col = LinearToGamma( vec4(col, 1.0), 0.8 ).rgb;//0.8

    return col;
}

#ifdef ENABLE_AUTO_VIEW
mat3 cameraAutoView( in sampler2D tex, out vec3 ro, out vec3 rd )
{
    float curTime = 10.0*time;
    ro = vec3(0.0, 0.0, -80.0-curTime);
    vec3 ta = vec3(0.0, 0.0, -90.0-curTime);
    ta = mix( ro + vec3(0.0, 1.0, 0.0), ta, smoothstep(1.0, 500.0, curTime) );
    ro.y = terrainL( tex, ro.xz ) + 30.0;//20.0 30.0
    ta.y = ro.y - 2.0;

    float fl = 1.2;//1.0;
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
    cameraAutoView( iChannel0, ro, rd );
#else
    @import ./ray;
#endif

    vec3 col = render( ro, rd );
    col = Vignetting( col, 0.5 );
    gl_FragColor = vec4( col, 1.0 );
}