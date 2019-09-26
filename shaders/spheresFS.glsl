@import ./include;

const float eps = 0.01;
const float OFFSET = eps * 2.0;//100.0;
const int STEPS = 64;

const vec3 lightDir = vec3( -0.48666426339228763, 0.8111071056538127, -0.3244428422615251 );

float random(vec2 p) { return fract(sin(mod(dot(p.xy, vec2(12.9898, 78.233)), 3.14))*43758.5453123); }

// distance functions
vec3 opRep( vec3 p, float interval )
{
    vec2 q = mod( p.xz, interval ) - interval * 0.5;
    return vec3( q.x, p.y, q.y );

    // vec3 r = vec3(interval);
    // return mod( p, r ) - 0.5 * r;
}

float sphereDist( vec3 p, float r )
{
    return length( opRep( p, 3.5*r ) ) - r;
}

float floorDist( vec3 p )
// floor.y = -1.0
{
    return dot(p, vec3( 0.0, 1.0, 0.0 ) ) + 1.0;
}

vec4 minVec4( vec4 a, vec4 b )
{
    return ( a.a < b.a ) ? a : b;
}

float checkerPattern( vec3 p )
{
    float u = 1.0 - floor( mod( p.x, 2.0 ) );
    float v = 1.0 - floor( mod( p.z, 2.0 ) );

    if ( ( u == 1.0 && v < 1.0 ) || ( u < 1.0 && v == 1.0 ) )
        return 0.5;
    else
        return 1.0;
}

vec3 hsv2rgb( vec3 c )
{
    vec4 K = vec4( 1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0 );
    vec3 p = abs( fract( c.xxx + K.xyz ) * 6.0 - K.www );
    return c.z * mix( K.xxx, clamp( p - K.xxx, 0.0, 1.0 ), c.y );
}

float sceneDist( vec3 p )
{
    float radius = 1.0;// + 0.2*sin(time*0.8);

    vec3 p1 = p;
    p1.y += cos(p.z + time*1.2) * sin(p.x + time*1.5);

    vec3 p2 = p;
    p2.y += cos(time) * sin(time);

    return min( sphereDist( p1, radius ), floorDist( p2 ) );
}

vec3 sceneColor( vec3 p )
{
    float radius = 1.0;
    vec3 sphereColor = hsv2rgb( vec3( (p.z + p.x)/9.0, 1.0, 1.0) );
    vec3 checkerColor = vec3(0.5) * checkerPattern( p );

    float d1 = sphereDist( p, radius );
    float d2 = floorDist( p );
    return mix( sphereColor, checkerColor, float(d1 > d2) );

    // return minVec4(
    //     // 3 * 6 / 2 = 9
    //     vec4( sphereColor, sphereDist( p, radius ) ),
    //     vec4( checkerColor, floorDist( p ) )
    // );
}

// vec3 getNormal( vec3 p )
// {
//     return normalize(vec3(
//         sceneDist(p + vec3( eps, 0.0, 0.0 ) ) - sceneDist(p + vec3( -eps, 0.0, 0.0 ) ),
//         sceneDist(p + vec3( 0.0, eps, 0.0 ) ) - sceneDist(p + vec3( 0.0, -eps, 0.0 ) ),
//         sceneDist(p + vec3( 0.0, 0.0, eps ) ) - sceneDist(p + vec3( 0.0, 0.0, -eps ) )
//     ));
// }
vec3 getNormal( vec3 pos, float eps )
{
    const vec3 v1 = vec3( 1.0,-1.0,-1.0);
    const vec3 v2 = vec3(-1.0,-1.0, 1.0);
    const vec3 v3 = vec3(-1.0, 1.0,-1.0);
    const vec3 v4 = vec3( 1.0, 1.0, 1.0);
    return normalize( v1 * sceneDist( pos + v1*eps ) +
        v2 * sceneDist( pos + v2*eps ) +
        v3 * sceneDist( pos + v3*eps ) +
        v4 * sceneDist( pos + v4*eps ) );
}
vec3 getNormal(vec3 pos)
{
    return getNormal(pos, 0.002);
}

float getShadow( vec3 ro, vec3 rd )
{
    float h = 0.0;
    float c = 0.0;
    float r = 1.0;
    float shadowCoef = 0.5;

    for ( float t = 0.0; t < 50.0; t++ )
    {
        h = sceneDist( ro + rd * c );
        if ( h < eps ) return shadowCoef;

        r = min( r, h * 16.0 / c );
        c += h;
    }

    return 1.0 - shadowCoef + r * shadowCoef;
}

vec3 getRayColor( vec3 ro, vec3 rd, out vec3 pos, out vec3 normal, out bool hit )
{
    // marching loop
    float dist;
    float depth = 0.0;
    pos = ro;

    for ( int i = 0; i < STEPS; i++ )
    {
        dist = sceneDist( pos );
        depth += dist;
        pos = ro + depth * rd;
        if ( abs(dist) < eps ) break;
    }

    // hit check and calc color
    vec3 color = vec3( 0.0 );
    if ( abs(dist) < eps )
    {
        normal = getNormal( pos );
        float diffuse = clamp( dot( lightDir, normal ), 0.1, 1.0 );
        float specular = pow( clamp( dot( reflect( lightDir, normal ), rd ), 0.0, 1.0 ), 10.0 );
        float shadow = getShadow( pos + normal * OFFSET, lightDir );
        color = ( sceneColor( pos ) * diffuse + vec3( 0.8 ) * specular ) * max( 0.5, shadow );
        hit = true;
    }

    // attenuation
    //return color - pow( clamp( 0.05 * depth, 0.0, 0.6 ), 2.0 );
    return color;
}

void main()
{
    @import ./ray;

    vec3 color = vec3( 0.0 );
    vec3 pos, normal;
    bool hit;
    float alpha = 1.0;

    for( int i = 0; i < 3; i++ )
    {
        color += alpha * getRayColor( ro, rd, pos, normal, hit );

        // compute fragDepth at the 1st hit...
        if( i == 0 && hit ) gl_FragDepth = getFragDepth( pos );

        // create a new ray for another hit...
        alpha *= 0.3;
        rd = normalize( reflect( rd, normal ) );
        ro = pos + normal * OFFSET; // if OFFSET = 0.0, ray does not start...

        if( !hit ) {
            gl_FragDepth = 1.0; // fragment will be killed
            break;
        }
    }

    gl_FragColor = vec4( color, 1.0 );
}