uniform sampler2D t_particle;
uniform float curTime;

void main()
{
  //vec4 col = vec4( 1.0, 0.3, 0.2, 1.0 );    // red flames
  //vec4 col = vec4( 0.07, 0.09, 0.52, 1.0 );   // blue flames
  vec4 col = vec4( 1.0 );
  vec4 tex = texture( t_particle, gl_PointCoord );
  gl_FragColor = vec4( col.rgb*tex.a, col.a*tex.a );
}


// vec4 getFlameColor(vec2 fragCoord)
// {
// 	float shift = 0.02;
// 	vec2 speed = vec2(1.5);

// 	const vec3 c1 = vec3(0.0470588235, 0.0471698113, 0.278431373);
// 	const vec3 c2 = vec3(0.854901961, 0.862745098, 0.0150993377);
// 	const vec3 c3 = vec3(0.2, 0.0, 0.0);
// 	const vec3 c4 = vec3(0.635294118, 0.00392156863, 0.00460732984);
// 	const vec3 c5 = vec3(3.1);
// 	const vec3 c6 = vec3(1.151);

//   // screen size
//   // large monitor(4K): 2775 1500
//   // small monitor(2K): 2066 1101
//   vec2 screenSize = vec2(2066.0, 1101.0);

// 	vec2 p = fragCoord.xy * 8.0 / screenSize.xx;
// 	float q = fbm(p - curTime * 0.1);
// 	vec2 r = vec2(fbm(p + q + curTime * speed.x - p.x - p.y), fbm(p + q - curTime * speed.y));
// 	vec3 c = mix(c1, c2, fbm(p + r)) + mix(c3, c4, r.x) - mix(c5, c6, r.y);
// 	float grad = fragCoord.y / screenSize.y;
// 	vec4 color = vec4(c * cos(shift * fragCoord.y / screenSize.y), 1.0);
// 	color.xyz *= 1.0-grad;
// 	return color;
// }
// void main()
// {
//   vec4 texColor = texture( t_particle, gl_PointCoord );//gl_PointCoord
// 	if( texColor.a < 0.5 ) discard;
//   //texColor.r += 0.5;

// 	vec3 fragColor = vec3(0.75, 0.5, 0.5);
//   vec4 flameColor = getFlameColor(gl_FragCoord.xy);//gl_FragCoord.xy
//   fragColor = mix(fragColor.xyz, flameColor.xyz, 0.8);

//   gl_FragColor = vec4( fragColor, 1.0 );
// }