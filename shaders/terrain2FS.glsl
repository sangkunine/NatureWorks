@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;
#define USE_WATER_1
//#define USE_CLOUDS_1
@import ./scene;

// textureMaps[0] = 'images/terrain/rockyTerrainDirt.jpg' (512 x 512)
// textureMaps[1] = 'images/raymarch/organic1.jpg' (1024 x 1024)
// textureMaps[2] = 'images/raymarch/grayNoise256.png' (256 x 256)
// textureMaps[3] = 'images/raymarch/lichen.jpg' (1024 x 1024)
#define iChannel0   textureMaps[0]
#define iChannel1   textureMaps[1]
#define iChannel2   textureMaps[2]
#define iChannel3   textureMaps[3]

// Choose...
#define ENABLE_AUTO_VIEW

const float FAR = 300.0;//300.0;
const vec3 SUN_LIGHT = normalize(vec3(-0.8, 0.4, -0.3));
const vec3 SUN_COLOR = vec3(1.0, 0.9, 0.85);
const float SKY_HEIGHT = 100.0;//300.0;
const vec3 WATER_COLOR = vec3(0.3, 0.4, 0.45);

// terrain land scales...
const float MYSTERY_LAND_SCALE = 1.65;
const float GRASS_LAND_SCALE = 1.75;
const float BUMPY_LAND_SCALE = 1.2;
const float STONE_LAND_SCALE = 1.9;

float skyMap( in vec3 p )
{
    return SKY_HEIGHT - p.y;
}

float grassLand( vec3 p )
{
    float h = triSmooth12( vec2(0.19, 0.31) + triSmooth12( p.xz/16.0 ) );
    //float h = triSmooth12( vec2(0.11, 0.99) + triSmooth12( p.xz/16.0 ) );
    //float h = triSmooth12( vec2(0.22, 0.88) + triSmooth12( p.xz/16.0 ) );
    h += triNoise13( p * 0.1 );
    return p.y - h * GRASS_LAND_SCALE;
}
float bumpyLand( vec3 p )
{
    float h = triNoise13( p * 0.1 );
    return p.y - h * BUMPY_LAND_SCALE;
}
float stoneLand( vec3 p )
{
    float h = triNoise13( p * 0.01 )*11.0; // 0.01 0.02 (for more complex stone)
    return p.y - h * STONE_LAND_SCALE;
}
float mysteryLand( vec3 p )
{
    float h = 0.0;
    h += triSmooth12( p.xz/16.0 )*0.66 + triSmooth12( p.xz/8.0 )*0.34;
    h += tri12( p.xz/2.0 )*0.23; p.xz = rotM2i * p.xz;
    h += tri12( p.xz*1.0 )*0.11; //p.xz = rotM2i * p.xz;
    h += triNoiseSmooth13( p * 0.01 )*11.5;
    return p.y - h * MYSTERY_LAND_SCALE;
}

vec2 terrainMap( in vec3 p )
{
    vec2 res = vec2( stoneLand(p), MATERIAL_TEXTURE0 );
    res = dUnion( res, vec2( bumpyLand(sTranslate(p, vec3(0.0, 4.0, 0.0))), MATERIAL_TEXTURE1 ) );
    return res;
}

float terrainL( in sampler2D tex, in vec3 p )
{
    return p.y - terrainMap( p ).x;
}

vec2 sceneMap( in vec3 p )
{
    vec2 res = vec2( skyMap(p), MATERIAL_SKY );
    res = dUnion( res, terrainMap(p) );
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

    vec3 col = vec3(0.0);

    if( tm.y == MATERIAL_SKY )
    {
        gl_FragDepth = 0.99;

        col = skyColor( SUN_LIGHT, rd, 0.0 );
        applyClouds( col, iChannel2, ro, rd );
	}

    else if( tm.y == MATERIAL_TEXTURE0 )
    {
        vec3 p = ro + rd * tm.x;
        gl_FragDepth = getFragDepth( p );

        // bump normal
        vec3 n = sceneNormal(p);
        //n = getBumpNormal(iChannel0, p*0.1, n, 0.1-0.05*tm.x/tmax );//0.075
        n = getBumpNormal(iChannel0, p*0.1, n, 0.075 );//0.075

        // cook shading
        vec3 diffuse, specular;
        float fre = pow( saturate(1.0 + dot(n, rd)), 5.0 );
        vec3 albedo = texCube( iChannel0, p, n ).xyz;
        albedo = mix(albedo, vec3(0.2,0.3,0.0), noise(iChannel2, p.xz*0.25));
        getCookShading( albedo, 0.01, 1.0, SUN_LIGHT, SUN_COLOR*3.0, n, rd, diffuse, specular );
        col = mix(diffuse, fre*vec3(1.0), saturate(n.y)) + specular;

        // curvature
        float curv = cheapCurvature(p)*0.9 + 0.1;
        col *= smoothstep(0.0, 0.8, curv);
    }
    else if( tm.y == MATERIAL_TEXTURE1 )
    {
        vec3 p = ro + rd * tm.x;
        gl_FragDepth = getFragDepth( p );

        // bump normal
        vec3 n = sceneNormal(p);
        //n = getBumpNormal(iChannel1, p*0.1, n, 0.1-0.05*tm.x/tmax );//0.075
        n = getBumpNormal(iChannel1, p*0.1, n, 0.075 );//0.075

        // cook shading
        vec3 diffuse, specular;
        float fre = pow( saturate(1.0 + dot(n, rd)), 5.0 );
        vec3 albedo = texCube( iChannel1, p, n ).xyz;
        albedo = mix(albedo, vec3(0.85, 0.4, 0.2), noise(iChannel2, p.xz*0.05));
        getCookShading( albedo, 0.01, 1.0, SUN_LIGHT, SUN_COLOR*3.0, n, rd, diffuse, specular );
        col = mix(diffuse, fre*vec3(0.7), saturate(n.y)) + specular;

        // curvature
        float curv = cheapCurvature(p)*0.9 + 0.1;
        col *= smoothstep(0.0, 0.8, curv);
    }

    // water surface...
    if( rd.y < 0.0 )
    {
        float waterHeight = 5.0;

        vec4 waterCol = waterColor( SUN_COLOR, SUN_LIGHT, WATER_COLOR, ro - vec3(0.0, waterHeight, 0.0), rd, vec4(0.0) );//cloudy = 0.0
        if( 0.0 < waterCol.w && waterCol.w < tm.x )
        {
            // foam on water...
            float t = (waterHeight-ro.y)/rd.y;
            vec2 uv = (ro + rd * t).xz;
            #if 0
                float sur = texture( iChannel3, 0.007*time * 0.06*uv ).x;
            #else
                float sur = texture( iChannel3, 0.06*uv ).x;
            #endif
            sur = smoothstep( 0.5, 1.0, sur )*0.5 + 0.5*sur*sur*smoothstep(0.2, 1.0, texture( iChannel2, 1.0*uv ).x);
            waterCol.rgb = mix( waterCol.rgb, vec3(2.5), 0.5*sur ); // foamCol = vec3(2.5)

            // sun specular...
            // float sunAmount = saturate( dot(SUN_LIGHT, reflect( rd, vec3(0.0,1.0,0.0) ) ) );
            // waterCol.rgb += 0.2*vec3(1.0,0.95,0.9)*pow(sunAmount,16.0);
            // waterCol.rgb += 0.5*vec3(1.0,0.95,0.9)*pow(sunAmount,96.0);

            col = mix(col, waterCol.rgb, saturate(1.1+rd.y));
        }
    }

    if( tm.y != MATERIAL_SKY )
    {
        // height-based fog density
        applyFog( col, SUN_COLOR, SUN_LIGHT, ro, rd, tm.x*1.4 );//1.2
        // constant fog density
        applyFog( col, SUN_COLOR, SUN_LIGHT, rd, 0.005, tm.x );//0.003
    }

    // post-processing
    col = FilmicToneMapping( col );
    col = col * 1.05 - 0.05;
    //col = col * 1.1 - 0.1;
    //col = LinearToGamma( vec4(col, 1.0), 0.8 ).rgb;//0.8

	return col;
}

#ifdef ENABLE_AUTO_VIEW
vec3 cpath( float t )
{
	vec3 p = vec3( 0.0, 0.0, 80.0 + t );//95.0
	float a = smoothstep(5.0, 20.0, t);
	p.xz += a*150.0 * cos( vec2(5.0,6.0) + 1.0*0.01*t );
	p.xz -= a*150.0 * cos( vec2(5.0,6.0) );
	p.xz += a* 50.0 * cos( vec2(0.0,3.5) + 6.0*0.01*t );
	p.xz -= a* 50.0 * cos( vec2(0.0,3.5) );
	return -p;
}

mat3 cameraAutoView( in sampler2D tex, out vec3 ro, out vec3 rd )
{
    float curTime = 2.0*time;//10.0
#if 0
    ro = vec3(0.0, 0.0, -80.0-curTime);
    vec3 ta = vec3(0.0, 0.0, -90.0-curTime);
#else
    ro = cpath( curTime );
	vec3 ta = cpath( 10.0 + curTime );
#endif
    ta = mix( ro + vec3(0.0, 1.0, 0.0), ta, smoothstep(1.0, 100.0, curTime) );
    ro.y = terrainL( tex, ro.xyz ) + 20.0;//30.0 20.0
    ta.y = ro.y - 5.0;

    float fl = 1.5;
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