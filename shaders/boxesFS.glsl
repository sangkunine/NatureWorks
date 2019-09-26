@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;
@import ./scene;

#define texChannel0 textureMaps[0]

//==============================================================================

// Choose...
#define ENABLE_AUTO_VIEW
#define ENABLE_SCENE_REPEAT
#define ENABLE_REFLECT_PASS

// Choose...
#define ENABLE_PHONG_SHADING
//#define ENABLE_COOK_SHADING

const vec3 lightCol = vec3(1.0, 0.97, 0.92) * 5.0;

vec2 sceneMap( vec3 p )
{
#ifdef ENABLE_SCENE_REPEAT
    p = mod(p, 1.0) - 0.5;
#endif

    //float d0 = fRoundBox(p, vec3(0.225,0.225,0.225), 0.1);// 0.1=radius (0.0=box, 0.5=sphere)
    //float d0 = fSphere(p, 0.25);
    float d0 = fBox(p, vec3(0.22, 0.22, 0.22));

    // float displace = 0.015*texture(texChannel0,(p.xy*p.xz)*2.0).x;
    // d0 -= displace;

    d0 = d0 - 0.2*fSinusoidalPlasma(p*2.0) + 0.1*fSinusoidalPlasma(p*4.0) - 0.05*fSinusoidalPlasma(p*8.0);
    
    return vec2(d0, 0.0);
}

vec3 getLighting( vec3 p, vec3 n, vec3 ld, vec3 rd, int reflectPass )
{
    // base color
    vec3 col = vec3(1.0);
    #ifdef ENABLE_SCENE_REPEAT
        vec3 voxPos = mod( p*0.5, 1.0 );
        vec3 spectumColor = vec3( sin(voxPos.x*7.0), 0.8*sin(voxPos.y*8.), 0.4*sin(voxPos.z*4.0) );
        if( (voxPos.x < 0.5) && (voxPos.z < 0.5) ) col = vec3(1.0, 0.85, 0.5) * (0.85 + 0.3*spectumColor);
        else if( (voxPos.x >= 0.5) && (voxPos.z >= 0.5) ) col = vec3(1.0, 0.85, 0.5) * (0.85 + 0.3*spectumColor);
    #endif

    // texture color
    col *= texCube(texChannel0, p, n, 4.0).xyz;

    // bump normal
    if( reflectPass == 0 ) n = getBumpNormal(texChannel0, p, n, 0.0075);

    // shadow & ao
    float shadow = (reflectPass == 0)? sceneShadow( p, ld ) : 1.0;
    float ao = sceneAO(p, n);

    // shading...
    #if defined (ENABLE_PHONG_SHADING)
        float shininess = 64.0;
        return getPhongShading( col, shininess, ld, lightCol, p, n, rd ) * (shadow * ao);

    #elif defined (ENABLE_COOK_SHADING)
        float metallic = 0.1;
        float roughness = 0.1;
        return getCookShading( col, metallic, roughness, ld, lightCol*2.0, n, rd ) * (shadow * ao);

    #endif
}

vec4 getSceneColor( vec3 ro, vec3 rd, float tmin, float tmax )
{
    //----------------------------------------------------
    // step 1: primary scene color: direct lighting pass
    //----------------------------------------------------

    gl_FragDepth = 0.99;

    // intersection
	float t = rayTracing( ro, rd, tmin, tmax, 0.75 );
	if( t > tmax ) return vec4(vec3(0.0), 1.0);

    // hitted
	vec3 p = ro + rd * t;
    vec3 n = sceneNormal( p, 0.2 );

    // light
    vec3 lp = vec3(ro.x, ro.y + 1.0, ro.z);
    vec3 ld = normalize(lp - p);

	vec3 sceneColor = getLighting( p, n, ld, rd, 0 ); // 0 = reflectPass ==> shadow will be computed
    gl_FragDepth = getFragDepth( p );

    //----------------------------------------------------
    // step 2: secondary scene color: reflection pass
    //----------------------------------------------------
    #ifdef ENABLE_REFLECT_PASS
        vec3 r = reflect(rd, n);
        t = rayTracing(p, r, 0.005, tmax, 0.75 ); // 0.005 (recommended: distance per one pixel)
        if( t < tmax ) {
            float refCoef = 0.35;
            p += r * t;
            n = sceneNormal( p, 0.2 );
            sceneColor += getLighting( p, n, ld, r, 1 ) * refCoef; // 1 = reflectPass ==> shadow will be computed
        }
    #endif

	return vec4( saturate(sceneColor), 1.0 );
}

mat2 rotMatrix( float angle )
{
	float c = cos( angle );
	float s = sin( angle );	
	return mat2( c, s, -s, c );
}
void setCamera( out vec3 ro, out vec3 rd )
{
    vec2 aspect = vec2(resolution.x/ resolution.y, 1.0);
	vec2 screenPos = (2.0*gl_FragCoord.xy / resolution.xy - 1.0) * aspect;

	float rad = 0.5;
	vec3 lookAt = vec3( 0.5,                       5.8*rad*cos(time*0.125), time       );
	vec3 camPos = vec3( 0.5+1.4*rad*sin(time*0.5), 5.8*rad*cos(time*0.125), -1.0 + time);

    vec3 forward = normalize(lookAt - camPos);
    vec3 right = normalize(vec3(forward.z, 0.0, -forward.x));
    vec3 up = normalize(cross(forward, right));
    float FOV = 0.5;
    ro = camPos;
    rd = normalize(forward + FOV*screenPos.x*right + FOV*screenPos.y*up);

    //rd.xz *= rotMatrix( PI*sin(-time*0.125)/2.0  ); // rotY ==> yawing sinusoidally (not frequently)
	rd.yz *= rotMatrix( PI*sin(-time*0.125)/6.0 ); // rotX ==> pitching (like nodding) sinusoidally
	rd.xy *= rotMatrix( PI*sin(-time*0.5)/4.0 );   // rotZ ==> rolling (like car wipers) sinusoidally
}

void main()
{
#ifdef ENABLE_AUTO_VIEW
    vec3 ro, rd;
    setCamera( ro, rd );
#else
    @import ./ray;
#endif

	float NEAR = 0.0;
	float FAR = 16.0;
    gl_FragColor = getSceneColor( ro, rd, NEAR, FAR );
}