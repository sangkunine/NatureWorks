@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;
#define USE_DUSTWIND_1
#define USE_LENSFLARES_1
#define USE_TERRAIN_3
@import ./scene;

// textureMaps[0] = heightTex0 ('images/raymarch/organic2.png')
// textureMaps[1] = heightTex1 ('images/raymarch/abstract1.jpg')
// textureMaps[2] = detailTex  ('images/raymarch/rgbaNoise256.png')

#define ENABLE_AUTO_VIEW
#define ENABLE_DUST_WIND
// #define ENABLE_SNOW
#define ENABLE_LENS_FLARES

const float FAR = 150.0;
const vec3 SUN_LIGHT = normalize(vec3(-1.0,0.19,0.4));
//const vec3 SUN_COLOR = vec3(1.0, 0.9, 0.85);
const vec3 SUN_COLOR = vec3(1.0, 0.9, 0.05);
const float SKY_HEIGHT = 500.0;//1000.0;

// vec3 SUN_COLOUR = vec3(1.1, 0.95, 0.85);
// vec3 FOG_COLOUR = vec3(0.48, 0.49, 0.53);
// vec3 SKY_BLUE_COLOUR = vec3(0.218, 0.325, 0.455);

float skyMap( in vec3 p )
{
    return SKY_HEIGHT - p.y;
}

vec2 sceneMap( in vec3 p )
{
	vec2 res = vec2( skyMap(p), MATERIAL_SKY );
	res = dUnion( res, vec2( terrainMap(textureMaps[0], textureMaps[1], textureMaps[2], p), MATERIAL_TERRAIN ) );
	return res;
}

void rayMinMax( in vec3 ro, in vec3 rd, out float tmin, out float tmax )
{
	tmin = 0.1;//1.0;
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
        
    vec2 tm = rayMarching( ro, rd, tmin, tmax );
	if( tm.x > tmax ) tm.y = MATERIAL_SKY;

	vec3 col = vec3(0.0);

	if( tm.y == MATERIAL_SKY )
	{
		gl_FragDepth = 0.99;
		col = skyColor( SUN_LIGHT, rd );
		applyClouds( col, ro, rd );
	}
	else if( tm.y == MATERIAL_TERRAIN )
    {
        vec3 p = ro + rd * tm.x;
		gl_FragDepth = getFragDepth( p );
		col = terrainColor( textureMaps[0], textureMaps[1], textureMaps[2], SUN_LIGHT, rd, p, tm.x );

		// col *= 1.0 when the sun is in front ===> as it is
		// col *= 0.5 when the sun is on the side ==> it gets dark
		col *= 1.0 - 0.5 * pow( 1.0-saturate(dot(rd,SUN_LIGHT)), 3.0);//0.25

		applyFog( col, SUN_COLOR, SUN_LIGHT, ro, rd, tm.x*1.2 );
		applyFog( col, SUN_COLOR, SUN_LIGHT, rd, 0.003, tm.x );//0.005

	#ifdef ENABLE_DUST_WIND
		float dustAmount = 0.2;//0.1;//0.01;
		float dustHeight = 10.0;//3.0;
		float windTurbulency = 0.3;
		applyDustWind( col, ro, rd, tm.x, dustAmount, dustHeight, windTurbulency );
	#endif
	}

    // sun glare
	col += sunScatter( SUN_LIGHT, rd );
	col = saturate(col);

	// color grading
#if 0
	col *= vec3(1.1, 1.0, 1.0); // almost linear
	col = col*col*(3.0-2.0*col); // col = smooth(0,1,col)
	col = pow( col, vec3(0.9,1.0,1.0) ); // red = sqr(red)

	col = mix( col, vec3(dot(col,vec3(0.333))), 0.4 ); // 0.6*col + 0.4*gray_of_col ==> gray-tone added
	col = col*0.5+0.5*col*col*(3.0-2.0*col); // almost linear
#endif

	col = col * 1.05 - 0.05;
	//col = col*col*(3.0-2.0*col);
	col = FilmicToneMapping( col );
	col = LinearToGamma( vec4(col, 1.0), 1.6 ).rgb;

#ifdef ENABLE_SNOW
	col += snowColor(0.2, 0.5);
#endif

	return col;
}

#ifdef ENABLE_AUTO_VIEW
vec3 cpath( float t )
{
	vec3 p = vec3( 0.0, 0.0, 95.0 + t );
	float a = smoothstep(5.0, 20.0, t);
	p.xz += a*150.0 * cos( vec2(5.0,6.0) + 1.0*0.01*t );
	p.xz -= a*150.0 * cos( vec2(5.0,6.0) );
	p.xz += a* 50.0 * cos( vec2(0.0,3.5) + 6.0*0.01*t );
	p.xz -= a* 50.0 * cos( vec2(0.0,3.5) );
	return -p;
}

mat3 cameraAutoView( in sampler2D heightTex0, out vec3 ro, out vec3 rd )
{
	float curTime = 0.5*(time);//1.0
#if 0
    ro = vec3(0.0, 0.0, -95.0-curTime);
    vec3 ta = vec3(0.0, 0.0, -105.0-curTime);
#else
	ro = cpath( curTime );
	vec3 ta = cpath( 2.0 + curTime );//10.0
#endif
    ta = mix( ro + vec3(0.0, 1.0, 0.0), ta, smoothstep(1.0, 10.0, curTime) );
	ro.y = terrainL( heightTex0, ro.xz ) + 0.1;//5.0
	ta.y = ro.y - 1.5;//5.0

    float fl = 1.0;//1.2;
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
    mat3 cam = cameraAutoView( textureMaps[0], ro, rd );
#else
    @import ./ray;
#endif

    vec3 col = render( ro, rd );

#ifdef ENABLE_LENS_FLARES
	vec2 xy = (-resolution.xy + 2.0*gl_FragCoord.xy)/resolution.y;
	#ifdef ENABLE_AUTO_VIEW
		applyLensFlares( col, cam, xy, SUN_COLOR, SUN_LIGHT, 0.0 );
	#else
		applyLensFlares( col, ro, cameraWorldMatrix, cameraProjectionMatrixInverse, xy, SUN_COLOR, SUN_LIGHT, 0.0 );
	#endif
#endif

	col = Vignetting( col, 0.5 );
	gl_FragColor = vec4( col, 1.0 );
}