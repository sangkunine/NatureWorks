@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;
#define USE_SEA_1
@import ./scene;

//==============================================================================

const vec3 SUN_LIGHT = normalize( vec3(0.4, 0.4, 0.48) );

vec2 sceneMap( in vec3 p ){ return vec2(0.0); }

void main(void)
{
    @import ./ray;

    // mixing = sky + sea
    float FAR = 1000.0;
    vec3 sky = skyColor( SUN_LIGHT, rd, 1.0 );
    vec3 sea = seaColor( SUN_LIGHT, ro, rd, FAR );
    vec3 color = mix( sky, sea, pow( smoothstep(0.0, -0.05, rd.y), 0.3) );

    // post-processing
    gl_FragColor = LinearToGamma( vec4(color, 1.0), 2.2 );
}