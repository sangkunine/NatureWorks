@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;
#define USE_WATER_1
#define USE_CLOUDS_1
#define USE_RAIN_1
#define USE_FLASH_1
#define USE_LENSFLARES_1
@import ./scene;

// rgbaNoise256.png
#define iChannel0   textureMaps[0]

//==============================================================================

const vec3 SUN_LIGHT = normalize(vec3(0.35, 0.14, 0.3));
const vec3 SUN_COLOR = vec3(1.0, 0.7, 0.55);
const vec3 WATER_COLOR = vec3(0.3, 0.4, 0.45);

vec2 sceneMap( in vec3 p ){ return vec2(0.0); }

void main()
{
	@import ./ray;

    float curTime = time*0.5;

	vec4 cloudy;
	cloudy.x = cos(curTime* 0.25)*0.5 + 0.5;		// cloud amount
    cloudy.yzw = flashColor( cloudy.x, curTime );	// flash in cloudy weather

	vec3 col;
	if( rd.y > 0.0 )
		col = skyCloudsColor( SUN_COLOR, SUN_LIGHT, ro, rd, cloudy );
	else
		col = waterColor( SUN_COLOR, SUN_LIGHT, WATER_COLOR, ro, rd, cloudy ).rgb;

	vec2 xy = ndc * vec2(resolution.x/resolution.y, 1.0); // xy = [-1.77, 1.77] x [-1, 1]
	applyLensFlares( col, ro, cameraWorldMatrix, cameraProjectionMatrixInverse, xy, SUN_COLOR, SUN_LIGHT, cloudy.x );
	applyRain( col, iChannel0, xy, cloudy, curTime );

	col = (col*col*(3.0-2.0*col));
	col *= 0.55 + 0.45 * Vignetting( col, 0.8 );
	gl_FragColor = LinearToGamma( vec4(col, 1.0), 2.2 );
}
