#ifndef RAYMARCH_UTILS
#define RAYMARCH_UTILS

//==============================================================================
// Basic functions
//==============================================================================

#define TONE_MAPPING_EXPOSURE 0.4

vec3 FilmicToneMapping( vec3 color ) {
    color *= TONE_MAPPING_EXPOSURE;
    return saturate( (color*(2.51*color + 0.03)) / (color*(2.43*color + 0.59) + 0.14) );
}

float pow2( const in float x ) { return x*x; }
float pow3( const in float x ) { return x*x*x; }
float pow4( const in float x ) { float x2 = x*x; return x2*x2; }
float pow5( const in float x ) { float x2 = x*x; return x2*x2*x; }

float mapLinear( float x, float x0, float x1, float y0, float y1 ) {
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}

float average( const in vec3 color ) { return dot( color, vec3( 0.3333 ) ); }

vec2 smoothstepd( float a, float b, float x )
// return smoothstep and its derivative
{
	if( x < a ) return vec2( 0.0, 0.0 );
	if( x > b ) return vec2( 1.0, 0.0 );
    float ir = 1.0 / (b - a);
    x = (x - a)*ir;
    return vec2( x*x*(3.0 - 2.0*x), 6.0*x*(1.0 - x)*ir );
}

float rand( in float x ) {
    return fract(sin(dot(vec2(x+47.49,38.2467/(x+2.3)), vec2(12.9898, 78.233)))*(43758.5453));
}

float rand( in vec2 uv ) {
	const float a = 12.9898, b = 78.233, c = 43758.5453;
	float dt = dot( uv.xy, vec2( a,b ) ), sn = mod( dt, PI );
	return fract(sin(sn) * c);
}

vec3 Dithering( in vec3 color ) {
    float grid_position = rand( gl_FragCoord.xy );
    vec3 dither_shift_RGB = vec3( 0.25 / 255.0, -0.25 / 255.0, 0.25 / 255.0 );
    dither_shift_RGB = mix( 2.0 * dither_shift_RGB, -2.0 * dither_shift_RGB, grid_position );
    return color + dither_shift_RGB;
}

vec3 Vignetting( in vec3 color, in float strength ) {
    // strength(0.0~3.0) = no-effect(0.0), weakly(0.5), normal(1.0), max-effect(5.0)
	vec2 uv = gl_FragCoord.xy / resolution.xy;
	vec2 offset = (uv - 0.5) * sqrt(2.0);
	float dist = dot(offset, offset);
	float shade = mix( 1.0, 1.0 - strength, dist );	
	return color * shade;
}

// uv = [0,1] x [0,1]
// xy = [-1.7,1.7] x [-1,1] (1.77 = aspect)
//
// vec2 xy = (-1.0 + 2.0*uv) * vec2(resolution.x/resolution.y, 1.0);
// vec2 uv = 0.5 + 0.5*xy*vec2(resolution.y/resolution.x, 1.0);
//
// vec2 uv = gl_FragCoord.xy/resolution.xy;
// vec2 xy = (-resolution.xy + 2.0*gl_FragCoord.xy)/resolution.y;

//==============================================================================
// Filtering
//==============================================================================

float checkerGradBox( in vec2 p )
// box-filtered checker-board (2D)
{
    // filter kernel
    vec2 w = fwidth(p) + 0.001;

    // analytical integral (box filter)
    vec2 i = 2.0*(abs(fract((p-0.5*w)*0.5)-0.5)-abs(fract((p+0.5*w)*0.5)-0.5))/w;

    // xor pattern
    return 0.5 - 0.5*i.x*i.y;
}

float gridGradBox( in vec2 p )
// box-filtered grid-board (2D)
{
    const float N = 10.0; // grid ratio

	// filter kernel
    vec2 w = fwidth(p) + 0.001;

	// analytic (box) filtering
    vec2 a = p + 0.5*w;                        
    vec2 b = p - 0.5*w;           
    vec2 i = (floor(a)+min(fract(a)*N,1.0) - floor(b)-min(fract(b)*N,1.0))/(N*w);

    //pattern
    return (1.0-i.x)*(1.0-i.y);
}

vec3 checkerColor( in vec3 p, in float ntile, in float intensity )
// ntile = number of tile per length
{
    float f = checkerGradBox( ntile * p.xz );
    return (1.0 - intensity) + f * vec3(intensity);
}

vec3 gridColor( in vec3 p, in float ntile, in float intensity )
// ntile = number of tile per length
{
    float f = gridGradBox( ntile * p.xz );
    return (1.0 - intensity) + f * vec3(intensity);
}

//==============================================================================
// Color
//==============================================================================

// Color mixing...
// col = mix(col1, col2, smoothstep(t1, t2, t))
// t < t1 ==> col = col1
// t > t2 ==> col = col2
// t1 < t < t2 ==> col = (1-w)*col1 + w*col2 (w = smoothstep(t1, t2, t))

float euclideanModulo( float n, float m )
{
    return mod((mod(n,m) + m), m);
}

float hue2rgb( float p, float q, float t )
{
    if ( t < 0.0 ) t += 1.0;
    if ( t > 1.0 ) t -= 1.0;
    if ( t < 1.0 / 6.0 ) return p + ( q - p ) * 6.0 * t;
    if ( t < 1.0 / 2.0 ) return q;
    if ( t < 2.0 / 3.0 ) return p + ( q - p ) * 6.0 * ( 2.0 / 3.0 - t );
    return p;
}

vec3 hsl2rgb( vec3 hsl )
// hsl = vec3(h, s, l)
// h,s,l ranges are in 0.0 - 1.0
{
    vec3 rgb;
    float h = euclideanModulo( hsl.x, 1.0 );
    float s = clamp( hsl.y, 0.0, 1.0 );
    float l = clamp( hsl.z, 0.0, 1.0 );
    if ( s == 0.0 )
        rgb = vec3(l);
    else {
        float p = (l <= 0.5) ? l * ( 1.0 + s ) : l + s - ( l * s );
        float q = ( 2.0 * l ) - p;
        rgb.r = hue2rgb( q, p, h + 1.0 / 3.0 );
        rgb.g = hue2rgb( q, p, h );
        rgb.b = hue2rgb( q, p, h - 1.0 / 3.0 );
    }
    return rgb;
}

//==============================================================================
// Texture
//==============================================================================

vec4 texCube( in sampler2D tex, in vec3 p, in vec3 n, in float k )
// p = hitted point(x, y, z)
// n = unit normal vector of its tangent plane
// k = 4.0
{
    vec3 m = pow( abs(n), vec3(k) );
	vec4 tx = texture( tex, p.yz );
	vec4 ty = texture( tex, p.zx );
	vec4 tz = texture( tex, p.xy );
	return (tx*m.x + ty*m.y + tz*m.z) / (m.x + m.y + m.z);
}
vec4 texCube( in sampler2D tex, in vec3 p, in vec3 n, in float k, in vec3 dpdx, in vec3 dpdy )
// p = hitted point(x, y, z)
// n = unit normal vector of its tangent plane
// k = 4.0
// dpdx = dp/dx = partial deriv of p w.r.t. window x
// dpdy = dp/dy = partial deriv of p w.r.t. window y
{
    vec3 m = pow( abs(n), vec3(k) );
	vec4 tx = textureGrad( tex, p.yz, dpdx.yz, dpdy.yz );
	vec4 ty = textureGrad( tex, p.zx, dpdx.zx, dpdy.zx );
	vec4 tz = textureGrad( tex, p.xy, dpdx.xy, dpdy.xy );
	return (tx*m.x + ty*m.y + tz*m.z) / (m.x + m.y + m.z);
}
vec4 texCube( in sampler2D tex, in vec3 p, in vec3 n )
{
    return texCube( tex, p, n, 4.0 );
}

// http://www.iquilezles.org/www/articles/texture/texture.htm
vec4 texSmooth( in sampler2D tex, in vec2 res, in vec2 uv )
// get smooth texture interpolation
// res = texture resolution
// uv = texture coordinates
{
	uv = uv*res + 0.5;
	vec2 iuv = floor( uv );
	vec2 fuv = fract( uv );
	uv = iuv + fuv*fuv*(3.0-2.0*fuv);
	uv = (uv - 0.5)/res;
	return texture( tex, uv );
}
vec4 texSmooth( in sampler2D tex, in vec2 res, in vec2 uv, in vec2 duvdx, in vec2 duvdy )
// get smooth texture interpolation
// res = texture resolution
// uv = texture coordinates
// duvdx = d(uv)/dx, duvdy = d(uv)/dy
{
	uv = uv*res + 0.5;
	vec2 iuv = floor( uv );
	vec2 fuv = fract( uv );
	uv = iuv + fuv*fuv*(3.0-2.0*fuv);
	uv = (uv - 0.5)/res;
	return textureGrad( tex, uv, duvdx, duvdy );
}

float getGrey( vec3 p ){ return p.x*0.299 + p.y*0.587 + p.z*0.114; }

vec3 getBumpNormal( in sampler2D tex, in vec3 p, in vec3 n, in float bumpFactor )
// bumpFactor = 0.0075, 0.075
{
    const float eps = 0.001;
    float baseVal = getGrey(texCube(tex, p, n).rgb);
    vec3 grad = vec3( getGrey(texCube(tex, vec3(p.x+eps, p.y, p.z), n).rgb) - baseVal,
                      getGrey(texCube(tex, vec3(p.x, p.y+eps, p.z), n).rgb) - baseVal,
                      getGrey(texCube(tex, vec3(p.x, p.y, p.z+eps), n).rgb) - baseVal )/eps;
    grad -= n * dot( n, grad );
    //return normalize( n - grad*bumpFactor );
    return normalize( n + grad*bumpFactor );
}

#endif // RAYMARCH_UTILS