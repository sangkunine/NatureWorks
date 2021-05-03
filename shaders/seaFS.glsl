@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;
#define USE_SEA_1
@import ./scene;

//==============================================================================

// Choose...
#define ENABLE_AUTO_VIEW

const vec3 SUN_LIGHT = normalize( vec3(0.4, 0.4, 0.48) );

vec2 sceneMap( in vec3 p ){ return vec2(0.0); }

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

    // camera ray direction
    float fl = 1.0;//1.5 2.0
    vec2 xy = (-resolution.xy + 2.0*gl_FragCoord.xy)/resolution.y;
    mat3 cam = cameraMatrix( ro, ta ); // cam[0] = cu, cam[1] = cv, cam[2] = cw
    rd = normalize( xy.x*cam[0] + xy.y*cam[1] + fl*cam[2] );
    return cam;
}
#endif

void main(void)
{
#ifdef ENABLE_AUTO_VIEW
    vec3 ro, rd;
    cameraAutoView( ro, rd );
#else
    @import ./ray;
#endif

    // mixing = sky + sea
    float FAR = 1000.0;
    vec3 sky = skyColor( SUN_LIGHT, rd, 1.0 );
    vec3 sea = seaColor( SUN_LIGHT, ro, rd, FAR );
    vec3 color = mix( sky, sea, pow( smoothstep(0.0, -0.05, rd.y), 0.3) );

    // post-processing
    gl_FragColor = LinearToGamma( vec4(color, 1.0), 1.2 ); // 2.2 => 1.2
}