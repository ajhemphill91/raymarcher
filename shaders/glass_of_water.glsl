////////////////////////////////////////////////////////////////
//
//             HELPER FUNCTIONS/MACROS
//
////////////////////////////////////////////////////////////////

#define PI 3.14159265
#define TAU (2.0*PI)
#define PHI (sqrt(5.0)*0.5 + 0.5)

// Sign function that doesn't return 0
float sgn(float x) {
	return (x<0.0)?-1.0:1.0;
}

vec2 sgn(vec2 v) {
	return vec2((v.x<0.0)?-1.0:1.0, (v.y<0.0)?-1.0:1.0);
}

float square (float x) {
	return x*x;
}

vec2 square (vec2 x) {
	return x*x;
}

vec3 square (vec3 x) {
	return x*x;
}

float lengthSqr(vec3 x) {
	return dot(x, x);
}

// Maximum/minumum elements of a vector
float vmax(vec2 v) {
	return max(v.x, v.y);
}

float vmax(vec3 v) {
	return max(max(v.x, v.y), v.z);
}

float vmax(vec4 v) {
	return max(max(v.x, v.y), max(v.z, v.w));
}

float vmin(vec2 v) {
	return min(v.x, v.y);
}

float vmin(vec3 v) {
	return min(min(v.x, v.y), v.z);
}

float vmin(vec4 v) {
	return min(min(v.x, v.y), min(v.z, v.w));
}

//Higher power'ed norms
float length8(vec2 p)
{
    p = p*p; p = p*p; p = p*p;
	  return pow(p.x + p.y, 1.0/8.0);
}

float length8(vec3 p)
{
    p = p*p; p = p*p; p = p*p;
    return pow(p.x + p.y + p.z, 1.0/8.0);
}

float length4(vec2 p)
{
    p = p*p; p = p*p;
    return pow(p.x + p.y, 1.0/4.0);
}

float length4(vec3 p)
{
    p = p*p; p = p*p;
    return pow(p.x + p.y + p.z, 1.0/8.0);
}

//Basic rotation matrices
mat3 Rx(float theta)
{
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(1, 0, 0),
        vec3(0, c, -s),
        vec3(0, s, c)
    );
}

mat3 Ry(float theta)
{
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, 0, s),
        vec3(0, 1, 0),
        vec3(-s, 0, c)
    );
}

mat3 Rz(float theta)
{
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, -s, 0),
        vec3(s, c, 0),
        vec3(0, 0, 1)
    );
}

////////////////////////////////////////////////////////////////
//
//             PRIMITIVE DISTANCE FUNCTIONS
//
////////////////////////////////////////////////////////////////
//
// Conventions:
//
// Everything that is a distance function is called fSomething.
// The first argument is always a point in 2 or 3-space called <p>.
// Unless otherwise noted, (if the object has an intrinsic "up"
// side or direction) the y axis is "up" and the object is
// centered at the origin.
//
////////////////////////////////////////////////////////////////

float fSphere(vec3 p, float r) {
	return length(p) - r;
}

// Plane with normal n (n is normalized) at some distance from the origin
float fPlane(vec3 p, vec3 n, float distanceFromOrigin) {
	return dot(p, n) + distanceFromOrigin;
}

// Cheap Box: distance to corners is overestimated
float fBoxCheap(vec3 p, vec3 b) { //cheap box
	return vmax(abs(p) - b);
}

// Box: correct distance to corners
float fBox(vec3 p, vec3 b) {
	vec3 d = abs(p) - b;
	return length(max(d, vec3(0.0))) + vmax(min(d, vec3(0.0)));
}

// Same as above, but in two dimensions (an endless box)
float fBox2Cheap(vec2 p, vec2 b) {
	return vmax(abs(p)-b);
}

float fBox2(vec2 p, vec2 b) {
	vec2 d = abs(p) - b;
	return length(max(d, vec2(0.0))) + vmax(min(d, vec2(0.0)));
}

//Box with rounded edges
float roundBox(vec3 p, vec3 halfFaceWidth, float cornerRadius)
{
    return length(max(abs(p)-halfFaceWidth,0.0))-cornerRadius;
}

// Endless "corner"
float fCorner (vec2 p) {
	return length(max(p, vec2(0.0))) + vmax(min(p, vec2(0.0)));
}

// Blobby ball object. You've probably seen it somewhere. This is not a correct distance bound, beware.
float fBlob(vec3 p) {
	p = abs(p);
	if (p.x < max(p.y, p.z)) p = p.yzx;
	if (p.x < max(p.y, p.z)) p = p.yzx;
	float b = max(max(max(
		dot(p, normalize(vec3(1.0))),
		dot(p.xz, normalize(vec2(PHI+1.0, 1.0)))),
		dot(p.yx, normalize(vec2(1.0, PHI)))),
		dot(p.xz, normalize(vec2(1.0, PHI))));
	float l = length(p);
	return l - 1.5 - 0.2 * (1.5 / 2.0)* cos(min(sqrt(1.01 - b / l)*(PI / 0.25), PI));
}

// Cylinder standing upright on the xz plane
float fCylinder(vec3 p, float r, float height) {
	float d = length(p.xz) - r;
	d = max(d, abs(p.y) - height);
	return d;
}

// Capsule: A Cylinder with round caps on both sides
float fCapsule(vec3 p, float r, float c) {
	return mix(length(p.xz) - r, length(vec3(p.x, abs(p.y) - c, p.z)) - r, step(c, abs(p.y)));
}

// Distance to line segment between <a> and <b>, used for fCapsule() version 2below
float fLineSegment(vec3 p, vec3 a, vec3 b) {
	vec3 ab = b - a;
	float t = saturate(dot(p - a, ab) / dot(ab, ab));
	return length((ab*t + a) - p);
}


// Capsule version 2: between two end points <a> and <b> with radius r
float fCapsule(vec3 p, vec3 a, vec3 b, float r) {
	return fLineSegment(p, a, b) - r;
}

// Torus in the XZ-plane
float fTorus(vec3 p, float smallRadius, float largeRadius) {
	return length(vec2(length(p.xz) - largeRadius, p.y)) - smallRadius;
}

// A circle line. Can also be used to make a torus by subtracting the smaller radius of the torus.
float fCircle(vec3 p, float r) {
	float l = length(p.xz) - r;
	return length(vec2(p.y, l));
}

// A circular disc with no thickness (i.e. a cylinder with no height).
// Subtract some value to make a flat disc with rounded edge.
float fDisc(vec3 p, float r) {
	float l = length(p.xz) - r;
	return l < 0.0 ? abs(p.y) : length(vec2(p.y, l));
}

// Hexagonal prism, circumcircle variant
float fHexagonCircumcircle(vec3 p, vec2 h) {
	vec3 q = abs(p);
	return max(q.y - h.y, max(q.x*sqrt(3.0)*0.5 + q.z*0.5, q.z) - h.x);
	//this is mathematically equivalent to this line, but less efficient:
	//return max(q.y - h.y, max(dot(vec2(cos(PI/3), sin(PI/3)), q.zx), q.z) - h.x);
}

// Hexagonal prism, incircle variant
float fHexagonIncircle(vec3 p, vec2 h) {
	return fHexagonCircumcircle(p, vec2(h.x*sqrt(3.0)*0.5, h.y));
}

// Cone with correct distances to tip and base circle. Y is up, 0 is in the middle of the base.
float fCone(vec3 p, float radius, float height) {
	vec2 q = vec2(length(p.xz), p.y);
	vec2 tip = q - vec2(0.0, height);
	vec2 mantleDir = normalize(vec2(height, radius));
	float mantle = dot(tip, mantleDir);
	float d = max(mantle, -q.y);
	float projected = dot(tip, vec2(mantleDir.y, -mantleDir.x));

	// distance to tip
	if ((q.y > height) && (projected < 0.0)) {
		d = max(d, length(tip));
	}

	// distance to base ring
	if ((q.x > radius) && (projected > length(vec2(height, radius)))) {
		d = max(d, length(q - vec2(radius, 0.0)));
	}
	return d;
}

////////////////////////////////////////////////////////////////
//
//                DOMAIN MANIPULATION OPERATORS
//
////////////////////////////////////////////////////////////////
//
// Conventions:
//
// Everything that modifies the domain is named pSomething.
//
// Many operate only on a subset of the three dimensions. For those,
// you must choose the dimensions that you want manipulated
// by supplying e.g. <p.x> or <p.zx>
//
// <inout p> is always the first argument and modified in place.
//
// Many of the operators partition space into cells. An identifier
// or cell index is returned, if possible. This return value is
// intended to be optionally used e.g. as a random seed to change
// parameters of the distance functions inside the cells.
//
// Unless stated otherwise, for cell index 0, <p> is unchanged and cells
// are centered on the origin so objects don't have to be moved to fit.
//
//
////////////////////////////////////////////////////////////////

// Rotate around a coordinate axis (i.e. in a plane perpendicular to that axis) by angle <a>.
// Read like this: R(p.xz, a) rotates "x towards z".
// This is fast if <a> is a compile-time constant and slower (but still practical) if not.
void pR(inout vec2 p, float a) {
	p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

// Shortcut for 45-degrees rotation
void pR45(inout vec2 p) {
	p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}

// Repeat space along one axis. Use like this to repeat along the x axis:
// <float cell = pMod1(p.x,5);> - using the return value is optional.
float pMod1(inout float p, float size) {
	float halfsize = size*0.5;
	float c = floor((p + halfsize)/size);
	p = mod(p + halfsize, size) - halfsize;
	return c;
}

// Same, but mirror every second cell so they match at the boundaries
float pModMirror1(inout float p, float size) {
	float halfsize = size*0.5;
	float c = floor((p + halfsize)/size);
	p = mod(p + halfsize,size) - halfsize;
	p *= mod(c, 2.0)*2.0 - 1.0;
	return c;
}

// Repeat the domain only in positive direction. Everything in the negative half-space is unchanged.
float pModSingle1(inout float p, float size) {
	float halfsize = size*0.5;
	float c = floor((p + halfsize)/size);
	if (p >= 0.0)
		p = mod(p + halfsize, size) - halfsize;
	return c;
}

// Repeat only a few times: from indices <start> to <stop> (similar to above, but more flexible)
float pModInterval1(inout float p, float size, float start, float stop) {
	float halfsize = size*0.5;
	float c = floor((p + halfsize)/size);
	p = mod(p+halfsize, size) - halfsize;
	if (c > stop) { //yes, this might not be the best thing numerically.
		p += size*(c - stop);
		c = stop;
	}
	if (c <start) {
		p += size*(c - start);
		c = start;
	}
	return c;
}

// Repeat around the origin by a fixed angle.
// For easier use, num of repetitions is use to specify the angle.
float pModPolar(inout vec2 p, float repetitions) {
	float angle = 2.0*PI/repetitions;
	float a = atan(p.y, p.x) + angle/2.0;
	float r = length(p);
	float c = floor(a/angle);
	a = mod(a,angle) - angle/2.0;
	p = vec2(cos(a), sin(a))*r;
	// For an odd number of repetitions, fix cell index of the cell in -x direction
	// (cell index would be e.g. -5 and 5 in the two halves of the cell):
	if (abs(c) >= (repetitions/2.0)) c = abs(c);
	return c;
}

// Repeat in two dimensions
vec2 pMod2(inout vec2 p, vec2 size) {
	vec2 c = floor((p + size*0.5)/size);
	p = mod(p + size*0.5,size) - size*0.5;
	return c;
}

// Same, but mirror every second cell so all boundaries match
vec2 pModMirror2(inout vec2 p, vec2 size) {
	vec2 halfsize = size*0.5;
	vec2 c = floor((p + halfsize)/size);
	p = mod(p + halfsize, size) - halfsize;
	p *= mod(c,vec2(2.0))*2.0 - vec2(1.0);
	return c;
}

// Same, but mirror every second cell at the diagonal as well
vec2 pModGrid2(inout vec2 p, vec2 size) {
	vec2 c = floor((p + size*0.5)/size);
	p = mod(p + size*0.5, size) - size*0.5;
	p *= mod(c,vec2(2.0))*2.0 - vec2(1.0);
	p -= size/2.0;
	if (p.x > p.y) p.xy = p.yx;
	return floor(c/2.0);
}

// Repeat in three dimensions
vec3 pMod3(inout vec3 p, vec3 size) {
	vec3 c = floor((p + size*0.5)/size);
	p = mod(p + size*0.5, size) - size*0.5;
	return c;
}

// Mirror at an axis-aligned plane which is at a specified distance <dist> from the origin.
float pMirror (inout float p, float dist) {
	float s = sgn(p);
	p = abs(p)-dist;
	return s;
}

// Mirror in both dimensions and at the diagonal, yielding one eighth of the space.
// translate by dist before mirroring.
vec2 pMirrorOctant (inout vec2 p, vec2 dist) {
	vec2 s = sgn(p);
	pMirror(p.x, dist.x);
	pMirror(p.y, dist.y);
	if (p.y > p.x)
		p.xy = p.yx;
	return s;
}

// Reflect space at a plane
float pReflect(inout vec3 p, vec3 planeNormal, float offset) {
	float t = dot(p, planeNormal)+offset;
	if (t < 0.0) {
		p = p - (2.0*t)*planeNormal;
	}
	return sgn(t);
}

//Grid repetition based on c which is the grid x,y,z spacing as c.x, c.y, c.z
vec3 opRep(vec3 p, vec3 c)
{
    return mod(p, c)-0.5*c;
}

////////////////////////////////////////////////////////////////
//
//             OBJECT COMBINATION OPERATORS
//
////////////////////////////////////////////////////////////////
//
// We usually need the following boolean operators to combine two objects:
// Union: OR(a,b)
// Intersection: AND(a,b)
// Difference: AND(a,!b)
// (a and b being the distances to the objects).
//
// The trivial implementations are min(a,b) for union, max(a,b) for intersection
// and max(a,-b) for difference. To combine objects in more interesting ways to
// produce rounded edges, chamfers, stairs, etc. instead of plain sharp edges we
// can use combination operators. It is common to use some kind of "smooth minimum"
// instead of min(), but we don't like that because it does not preserve Lipschitz
// continuity in many cases.
//
// Naming convention: since they return a distance, they are called fOpSomething.
// The different flavours usually implement all the boolean operators above
// and are called fOpUnionRound, fOpIntersectionRound, etc.
//
// The basic idea: Assume the object surfaces intersect at a right angle. The two
// distances <a> and <b> constitute a new local two-dimensional coordinate system
// with the actual intersection as the origin. In this coordinate system, we can
// evaluate any 2D distance function we want in order to shape the edge.
//
// The operators below are just those that we found useful or interesting and should
// be seen as examples. There are infinitely more possible operators.
//
// They are designed to actually produce correct distances or distance bounds, unlike
// popular "smooth minimum" operators, on the condition that the gradients of the two
// SDFs are at right angles. When they are off by more than 30 degrees or so, the
// Lipschitz condition will no longer hold (i.e. you might get artifacts). The worst
// case is parallel surfaces that are close to each other.
//
// Most have a float argument <r> to specify the radius of the feature they represent.
// This should be much smaller than the object size.
//
// Some of them have checks like "if ((-a < r) && (-b < r))" that restrict
// their influence (and computation cost) to a certain area. You might
// want to lift that restriction or enforce it. We have left it as comments
// in some cases.
//
// usage example:
//
// float fTwoBoxes(vec3 p) {
//   float box0 = fBox(p, vec3(1));
//   float box1 = fBox(p-vec3(1), vec3(1));
//   return fOpUnionChamfer(box0, box1, 0.2);
// }
//
////////////////////////////////////////////////////////////////

//ONLY SDF a and b both are
vec2 intersectSDF(vec2 a, vec2 b)
{
    //max(a.x,b.x)
    return (a.x>b.x) ? a : b;
}

//ANYWHERE SDF a and b both are
vec2 unionSDF(vec2 a, vec2 b)
{
    //min(a.x,b.x)
    return (a.x<b.x) ? a : b;
}

float unionSDF(float a, float b)
{
		return min(a, b);
}

//Union with rounded corners
float smoothUnionSDF(float a, float b, float k)
{
    float res = exp(-k*a)+exp(-k*b);
    return -log(res)/k;
}

//Same idea, different implementation
float polynomialUnionSDF(float a, float b, float k)
{
    float h = clamp(0.5+0.5*(b-a)/k, 0.0, 1.0);
    return mix(b, a, h) - k*h*(1.0-h);
}

//Where SDF a IS AND b ISNT
vec2 diffSDF(vec2 a, vec2 b)
{
    return (a.x<-b.x) ? a : b;
}

float diffSDF(float a, float b)
{
		return max(a, -b);
}

// The "Chamfer" flavour makes a 45-degree chamfered edge (the diagonal of a square of size <r>):
float fOpUnionChamfer(float a, float b, float r) {
	return min(min(a, b), (a - r + b)*sqrt(0.5));
}

// Intersection has to deal with what is normally the inside of the resulting object
// when using union, which we normally don't care about too much. Thus, intersection
// implementations sometimes differ from union implementations.
float fOpIntersectionChamfer(float a, float b, float r) {
	return max(max(a, b), (a + r + b)*sqrt(0.5));
}

// Difference can be built from Intersection or Union:
float fOpDifferenceChamfer (float a, float b, float r) {
	return fOpIntersectionChamfer(a, -b, r);
}

// The "Round" variant uses a quarter-circle to join the two objects smoothly:
float fOpUnionRound(float a, float b, float r) {
	vec2 u = max(vec2(r - a,r - b), vec2(0.0));
	return max(r, min (a, b)) - length(u);
}

float fOpIntersectionRound(float a, float b, float r) {
	vec2 u = max(vec2(r + a,r + b), vec2(0.0));
	return min(-r, max (a, b)) + length(u);
}

float fOpDifferenceRound (float a, float b, float r) {
	return fOpIntersectionRound(a, -b, r);
}

// The "Columns" flavour makes n-1 circular columns at a 45 degree angle:
float fOpUnionColumns(float a, float b, float r, float n) {
	if ((a < r) && (b < r)) {
		vec2 p = vec2(a, b);
		float columnradius = r*sqrt(2.0)/((n-1.0)*2.0+sqrt(2.0));
		pR45(p);
		p.x -= sqrt(2.0)/2.0*r;
		p.x += columnradius*sqrt(2.0);
		if (mod(n,2.0) == 1.0) {
			p.y += columnradius;
		}
		// At this point, we have turned 45 degrees and moved at a point on the
		// diagonal that we want to place the columns on.
		// Now, repeat the domain along this direction and place a circle.
		pMod1(p.y, columnradius*2.0);
		float result = length(p) - columnradius;
		result = min(result, p.x);
		result = min(result, a);
		return min(result, b);
	} else {
		return min(a, b);
	}
}

float fOpDifferenceColumns(float a, float b, float r, float n) {
	a = -a;
	float m = min(a, b);
	//avoid the expensive computation where not needed (produces discontinuity though)
	if ((a < r) && (b < r)) {
		vec2 p = vec2(a, b);
		float columnradius = r*sqrt(2.0)/n/2.0;
		columnradius = r*sqrt(2.0)/((n-1.0)*2.0+sqrt(2.0));

		pR45(p);
		p.y += columnradius;
		p.x -= sqrt(2.0)/2.0*r;
		p.x += -columnradius*sqrt(2.0)/2.0;

		if (mod(n,2.0) == 1.0) {
			p.y += columnradius;
		}
		pMod1(p.y,columnradius*2.0);

		float result = -length(p) + columnradius;
		result = max(result, p.x);
		result = min(result, a);
		return -min(result, b);
	} else {
		return -m;
	}
}

float fOpIntersectionColumns(float a, float b, float r, float n) {
	return fOpDifferenceColumns(a,-b,r, n);
}

// The "Stairs" flavour produces n-1 steps of a staircase:
// much less stupid version by paniq
float fOpUnionStairs(float a, float b, float r, float n) {
	float s = r/n;
	float u = b-r;
	return min(min(a,b), 0.5 * (u + a + abs ((mod (u - a + s, 2.0 * s)) - s)));
}

// We can just call Union since stairs are symmetric.
float fOpIntersectionStairs(float a, float b, float r, float n) {
	return -fOpUnionStairs(-a, -b, r, n);
}

float fOpDifferenceStairs(float a, float b, float r, float n) {
	return -fOpUnionStairs(-a, b, r, n);
}

// Similar to fOpUnionRound, but more lipschitz-y at acute angles
// (and less so at 90 degrees). Useful when fudging around too much
// by MediaMolecule, from Alex Evans' siggraph slides
float fOpUnionSoft(float a, float b, float r) {
	float e = max(r - abs(a - b), 0.0);
	return min(a, b) - e*e*0.25/r;
}

// produces a cylindical pipe that runs along the intersection.
// No objects remain, only the pipe. This is not a boolean operator.
float fOpPipe(float a, float b, float r) {
	return length(vec2(a, b)) - r;
}

// first object gets a v-shaped engraving where it intersect the second
float fOpEngrave(float a, float b, float r) {
	return max(a, (a + r - abs(b))*sqrt(0.5));
}

// first object gets a capenter-style groove cut out
float fOpGroove(float a, float b, float ra, float rb) {
	return max(a, min(a + ra, rb - abs(b)));
}

// first object gets a capenter-style tongue attached
float fOpTongue(float a, float b, float ra, float rb) {
	return min(a, max(a - ra, abs(b) - rb));
}

///////////////////////////////////////////////////////////////////////////////
//
//													MISC USEFUL FUNCTIONS
//
///////////////////////////////////////////////////////////////////////////////

//SDFs of the returned value will be twisted around an axis
vec3 opTwist(vec3 p, float amt)
{
    float c = cos(amt*p.y);
    float s = sin(amt*p.y);
    mat2  m = mat2(c,-s,s,c);
    return vec3(m*p.xz, p.y);
}

//Same as twist but with a bendy action
vec3 opBend(vec3 p, float amt)
{
    float c = cos(amt*p.y);
    float s = sin(amt*p.y);
    mat2  m = mat2(c,-s,s,c);
    return vec3(m*p.xy, p.z);
}

//Example of displacement function, you can just add this to the SDF result
float sinDisplace(vec3 p)
{
    return (sin(iGlobalTime*p.x)+sin(iGlobalTime*p.y)+sin(iGlobalTime*p.z));
}

//Noise functions, from I. Quilez
vec2 random2(vec2 st)
{
    st = vec2( dot(st,vec2(127.1,311.7)),
              dot(st,vec2(269.5,183.3)) );
    return -1.0 + 2.0*fract(sin(st)*43758.5453123);
}

float noise(vec2 st)
{
    vec2 i = floor(st);
    vec2 f = fract(st);

    vec2 u = f*f*(3.0-2.0*f);

    return mix( mix( dot( random2(i + vec2(0.0,0.0) ), f - vec2(0.0,0.0) ),
                     dot( random2(i + vec2(1.0,0.0) ), f - vec2(1.0,0.0) ), u.x),
                mix( dot( random2(i + vec2(0.0,1.0) ), f - vec2(0.0,1.0) ),
                     dot( random2(i + vec2(1.0,1.0) ), f - vec2(1.0,1.0) ), u.x), u.y);
}

//Basic/naive simplex noise implementation, should be looked at again soon (ajh 6/2017)
vec3 simplexNoise()
{
    vec3 z = vec3(0.0);
    vec2 uv = gl_FragCoord.st/iResolution.xy;
    float r = 115.47*uv.x;
    vec2 s = fract(vec2(r, uv.y + r*0.5));
    if(s.x > s.y)
    {
      z.xy = 1.0-vec2(s.x, s.y-s.x);
      z.z = s.y;
    }
    else
    {
      z.yz = 1.0-vec2(s.x-s.y, s.y);
      z.x = s.x;
    }

    return fract(z);
}

//Oscillates at global time from lower bound to upper bound, both l and u should > 0
float sint(float l, float u)
{
    return (sin(iGlobalTime)/2.0+0.5)*u + l;
}

float cost(float l, float u)
{
    return (cos(iGlobalTime)/2.0+0.5)*u + l;
}

//Pack a float value as a vec3. Useful for distance to color(RGB)
vec3 unpackColor(float f)
{
    vec3 color;
    color.b = floor(f / 256.0 / 256.0);
    color.g = floor((f - color.b * 256.0 * 256.0) / 256.0);
    color.r = floor(f - color.b * 256.0 * 256.0 - color.g * 256.0);
    // now we have a vec3 with the 3 components in range [0..255]. Let's normalize it!
    return color / 255.0;
}

///////////////////////////////////////////////////////////////////////////////
//
//			END OF SIGNED DISTANCE FUNCTIONS, OPERATORS, UTILITIES
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
//															RAYTRACER
//
///////////////////////////////////////////////////////////////////////////////
//
//	map(p) defines the main scene via distance functions.
//		p is usually given as (eye + z*v), which gives
//				a point in space z units in the v direction from eye (the camera)
//		The return value is a vec2, x component is the result of SDFs on p,
//		and the y component can be assigned values as IDs, to allow for varying
//		shading methods, etc once the ray has been evaluated
//
///////////////////////////////////////////////////////////////////////////////

struct Ray_t
{
    vec3 ro;	//Ray origin
    vec3 rd;	//Ray direction
		float totalDistance;
    float fStartDistance;
    float fLength;
		bool wasTransmitted;	//Signifies this ray is to compute transmission
		bool active;					//Whether to skip this ray
};


struct Material_t
{
		vec3 diffuseColor;
		vec3 ambientColor;
		float reflections;
		float maxReflections;
    float transparency;	//Domain [0, 1], 1 is fully transparent, 0 = opaque
		float smoothness;	//Fresnel approximation co-effecient
    float n;	//Refractive index
};

struct Hit_t
{
		float distTravelled;	//Distance incident ray travelled to reach surface
    vec3 pos;							//Point of intersection
		vec3 N;								//Normal vector to intersected surface
		vec3 R;								//Reflection vector
		vec3 T;								//Transmission vector
		float materialID;
    //Material_t mat;				//Material at this point
};

struct Light_t
{
		vec3 ro;
		vec3 rd;
		vec3 color;
		float intensity;
		bool directional;
};

#define NUM_LIGHTS 1
Light_t _lights[NUM_LIGHTS];

struct Surface_t
{
    vec3 vNormal;
    vec3 cReflection;
    vec3 cTransmission;
};

struct C_Shading
{
    vec3 cDiffuse;
    vec3 cSpecular;
};

struct C_PointLight
{
    vec3 vPos;
    vec3 cColour;
};

struct C_DirectionalLight
{
    vec3 vDir;
    vec3 cColour;
};

#define EPSILON 								0.001
#define NEAR_CLIP_PLANE         1.0
#define FAR_CLIP_PLANE         	100.0
#define MAX_RAYCAST_STEPS       150
#define STEP_DAMPING            0.9
#define DIST_EPSILON            0.0001
#define MAX_RAY_BOUNCES         3.0
#define MAX_RAY_TRANSMISSIONS		3.0

#define MAX_SHADOW_DIST         10.0
#define SHADOW_K								32.0

#define AMBIENT_COLOR           vec3(0.00, 0.0, 0.00)
#define LIGHT_COLOR             vec3(0.8, 0.9, 0.8)
#define SPEC_COLOR              vec3(0.8, 0.9, 0.8)

#define SPEC_POWER             	64.0		//Incident light transmitted power, I

#define FOG_DENSITY             0.01

#define REFLECTIVITY						0.9

#define CAM_DIST                15.0
#define CAM_H                   1.0
#define CAM_FOV_FACTOR 					1.0
#define LOOK_AT_H               1.0
#define LOOK_AT 								vec3(0.0, LOOK_AT_H, 0.0)

#define MTL_BACKGROUND          -1.0
#define MTL_GLASS								1.0
#define MTL_TABLE								2.0
#define MTL_OJ									3.0
#define MTL_EMISSION						200.0

#define CLR_BACKGROUND 					vec3(0.1, 0.1, 0.15)

//#define ALLOW_POINT_LIGHTING

vec2 map(vec3 p)
{
    vec2 scene;

		p = Rz(PI/-4.)*p;
		p = Ry(iGlobalTime/4.)*p;
		vec3 pGlass = p + vec3(2.0, -1.2, 2.0);
		float table = fBox(p + vec3(0.0, 4.0, 0.0), vec3(15.0, 1.0, 15.0));
		float glass = fCylinder(pGlass, 3.0, 4.0);
		float cavity = fCylinder(pGlass - vec3(0.0, 0.35, 0.0), 2.9, 3.75);
		glass = fOpDifferenceRound(glass, cavity, 0.1);
		float oj = fCylinder(pGlass + vec3(0.0, 1.65, 0.0), 2.9, 2.0);

		scene = vec2(table, MTL_TABLE);
		scene = unionSDF(scene, vec2(glass, MTL_GLASS));
		scene = unionSDF(scene, vec2(oj, MTL_OJ));

    return scene;
}

vec3 normal(in vec3 p)
{
    vec2 d = vec2(EPSILON, 0.0);
    return normalize(vec3(
        map(p + d.xyy).x - map(p - d.xyy).x,
        map(p + d.yxy).x - map(p - d.yxy).x,
        map(p + d.yyx).x - map(p - d.yyx).x));
}

vec2 rayMarch(Ray_t ray)
{
    float t = NEAR_CLIP_PLANE;
    float materialID = MTL_BACKGROUND;
    for (int i = 0; i < MAX_RAYCAST_STEPS; i++)
		{
        vec2 d = map(ray.ro + ray.rd*t);
        if (d.x < DIST_EPSILON || t > FAR_CLIP_PLANE) break;
        t += d.x*STEP_DAMPING;
        materialID = d.y;
    }

    if (t > FAR_CLIP_PLANE) materialID = MTL_BACKGROUND;
    return vec2(t, materialID);
}

vec3 applyFog(vec3 col, float dist)
{
    return mix(col, CLR_BACKGROUND, 1.0 - exp(-FOG_DENSITY*dist*dist));
}

float shadow(vec3 ro, vec3 rd)
{
    float shadowAmt = 1.0;
    float t = NEAR_CLIP_PLANE;
    for (int i = 0; i < 16; i++)
		{
        float d = map(ro + rd*t).x*STEP_DAMPING;
				shadowAmt = min(shadowAmt, SHADOW_K*d/t);
        t += clamp(d, 0.02, 0.1);
        if (d < 0.001 || t > MAX_SHADOW_DIST) break;
    }

    return clamp(shadowAmt, 0.0, 1.0);
}

Material_t getMaterial(float matID)
{
		Material_t m;
		m.reflections = 0.0;
		m.maxReflections = 0.0;
		m.transparency = 0.0;
		m.smoothness = 0.0;
		m.n = 1.0;
		m.ambientColor = AMBIENT_COLOR;

		if(matID == MTL_GLASS)
		{
				m.diffuseColor = vec3(0.5);
				m.maxReflections = MAX_RAY_BOUNCES;
				m.smoothness = 0.9;
				m.transparency = 0.9;
				m.n = 1.1;
		}
		else if(matID == MTL_TABLE)
		{
				m.diffuseColor = vec3(0.2, 0.1, 0.1);
				m.smoothness = 0.3;
		}
		else if(matID == MTL_OJ)
		{
				m.diffuseColor = vec3(0.8, 0.5, 0.0);
				m.maxReflections = MAX_RAY_BOUNCES;
				m.smoothness = 1.0;
				m.transparency = 0.5;
				m.n = 1.3;
		}
		else if(matID == MTL_EMISSION)
		{
				m.diffuseColor = vec3(0.25);
				m.smoothness = 0.4;
				m.transparency = 0.9;
				m.n = 1.33;

		}
		else
		{
				m.diffuseColor = CLR_BACKGROUND;
				m.smoothness = 1.0;

		}

    return m;
}

float blinnPhongIntensity(vec3 incidentDir, vec3 lightDir, vec3 N, float smoothness)
{
    vec3 H = normalize(lightDir - incidentDir);		//Half-way vector
    float NdotH = max(0.0, dot(H, N));

    float specPower = exp2(4.0 + 6.0 * smoothness);
    float specIntensity = (specPower + 2.0) * 0.125;

    return pow(NdotH, specPower) * specIntensity;
}

vec3 computeColor(Ray_t ray, Hit_t hit, Material_t mat, Ray_t eye)
{
		//Evaluate the diffuse and specular colors
		vec3 res = vec3(0.0);
		for (int j = 0; j < NUM_LIGHTS; j++)
		{
				Ray_t lightRay;

				lightRay.rd = _lights[j].rd;
				if(!_lights[j].directional)
				{
						//Light direction from a point light
						Ray_t toLight;
						toLight.ro = hit.pos;
						toLight.rd = _lights[j].ro - hit.pos;

						//Check to be sure the light doesn't intersect another object on the way
						if(rayMarch(toLight).y != MTL_BACKGROUND)
						{
								continue;
						}

						lightRay.rd = normalize(_lights[j].ro - hit.pos);
						lightRay.ro = _lights[j].ro;

				}

				float lambertian = clamp(dot(lightRay.rd, hit.N), 0.0, 1.0);
				lambertian *=
					_lights[j].intensity / (ray.totalDistance * ray.totalDistance);

				lambertian *= shadow(hit.pos, lightRay.rd);

				float bpi = 0.0;
				//if(!(ray.wasTransmitted && !_lights[j].directional))
				{
					bpi = blinnPhongIntensity(eye.rd, lightRay.rd, hit.N, mat.smoothness) *
						_lights[j].intensity / (ray.totalDistance * ray.totalDistance);	//eye.rd is used here as the view vector
				}



				//vec3 fresColor = SPEC_COLOR + (1.0 - SPEC_COLOR)*pow((1.0 - specular), 5.0); //Schlick approximation for fresnel speculars
				lambertian *= shadow(ray.ro, lightRay.rd);

				vec3 col = mat.ambientColor + mat.diffuseColor*(_lights[j].color*lambertian + _lights[j].color*bpi);

				//col += applyFog(col, hit.distTravelled);
				res += col*(1.0 - mat.transparency); //Blend this ray's result to the pixel value, could also be colored by reflection/transmission
		}

		return res;
}

Hit_t computeHit(Ray_t ray, vec2 map)
{
		Hit_t hit;

		ray.totalDistance += map.x;
		hit.distTravelled = ray.totalDistance;
		hit.materialID = map.y;
		hit.pos = ray.ro + (ray.rd*map.x);
		hit.N = normal(hit.pos);
		hit.R = reflect(ray.rd, hit.N);


		//hit.T = refract(ray.rd, hit.N, 1.0/hit.mat.n);
		//TODO: The n1=1.0 is not right, should hold onto an index of refraction in the ray struct for the previous medium's index of refraction

		return hit;
}

vec3 computeTransmissions(Ray_t incident, Hit_t hit, Material_t mat)
{
		vec3 res = vec3(0.0);

		//Using the inverse square law to assume a reasonable cutoff for the recursion
		if(incident.totalDistance > 10000.0) return res;

		for (float i = 0.0; i < MAX_RAY_TRANSMISSIONS; i++)
		{
				if(mat.transparency == 0.0) break;

				//Compute the transmitted ray's effect
				Ray_t tx;
				tx.ro = hit.pos + hit.T*DIST_EPSILON;
				tx.rd = hit.T;

				for (float j = 0.0; j < MAX_RAY_BOUNCES; j++)
				{
						vec2 txResult = rayMarch(tx);
						Hit_t txHit = computeHit(tx, txResult);
						Material_t txMat = getMaterial(txHit.materialID);

						//Further transmission
						//res += computeTransmissions(tx, txHit);

						//Further reflection
						//res += computeColor(tx, txHit, txMat);

						//Break if the matrial shouldn't reflect or if the ray is glancing (large angle between normal and ray)
		        if (txMat.maxReflections == 0.0 || abs(dot(txHit.N, tx.rd)) > 0.99 ||
								txHit.materialID == MTL_BACKGROUND) break;	//TODO: should be dot > 0.99 or so

						//Otherwise set up the reflected ray and loop on
		        tx.ro = txHit.pos + txHit.R*DIST_EPSILON;
		        tx.rd = txHit.R;
		    }

		}

		return res;

}

vec3 render(Ray_t ray)
{
    vec3 resCol = vec3(0.0);
    float alpha = 1.0;

		Ray_t secondary = ray;

    for (float i = 0.0; i < MAX_RAY_BOUNCES; i++)
		{
        vec2 marchingResult = rayMarch(secondary);
				Hit_t hit = computeHit(secondary, marchingResult);
				Material_t mat = getMaterial(hit.materialID);

				//Transmitted component of light
				//resCol += computeTransmissions(secondary, hit);

				//Direct lighting and reflection components
				resCol += computeColor(secondary, hit, mat, ray);

				//Break if the matrial shouldn't reflect or if the ray is glancing (large angle between normal and ray)
        if (mat.maxReflections == 0.0 || abs(dot(hit.N, secondary.rd)) < 0.01 ||
						hit.materialID == MTL_BACKGROUND) break;	//TODO: should be dot > 0.99 or so

				//Otherwise set up the reflected ray and loop on
        secondary.ro = hit.pos + hit.R*DIST_EPSILON;
        secondary.rd = hit.R;
    }

		//Gamma correction
		//resCol = pow(resCol, vec3(1.0/2.2));

    return vec3(clamp(resCol, 0.0, 1.0));
}

vec3 getRayDir(vec3 camPos, vec3 viewDir, vec2 pixelPos)
{
    vec3 camRight = normalize(cross(viewDir, vec3(0.0, 1.0, 0.0)));
    vec3 camUp = normalize(cross(camRight, viewDir));
    return normalize(pixelPos.x*camRight + pixelPos.y*camUp + CAM_FOV_FACTOR*viewDir);
}

void main(void)
{
    vec2 q = gl_FragCoord.xy/iResolution.xy;
    vec2 p = -1.0+2.0*q;
    p.x *= iResolution.x/iResolution.y;

    float ang = 0.0;

    vec3 camPos = vec3(CAM_DIST*cos(ang), CAM_H, CAM_DIST*sin(ang));
		Ray_t ray;
		ray.totalDistance = 0.0;
		ray.ro = camPos;
    ray.rd = getRayDir(camPos, normalize(LOOK_AT - camPos), p);

		//Set up the lights
		_lights[0].rd = -ray.rd;	//Headlight on the camera
		_lights[0].directional = true;
		_lights[0].color = LIGHT_COLOR;
		_lights[0].intensity = 150.0;
/*
		_lights[1].ro = vec3(0.0, -10.0, 0.0);
		_lights[1].directional = false;		//Point light
		_lights[1].color = vec3(0.4, 0.0, 0.7);
		_lights[1].intensity = 10.0;

		_lights[2].ro = vec3(15.0, -10.0, 0.0);
		_lights[2].directional = false;
		_lights[2].color = vec3(0.9, 0.4, 0.0);
		_lights[2].intensity = 20.0;

		_lights[3].ro = vec3(0.0, -10.0, 15.0);
		_lights[3].directional = false;
		_lights[3].color = vec3(0.1, 0.7, 0.3);
		_lights[3].intensity = 10.0;

		_lights[4].ro = vec3(15.0, -10.0, 15.0);
		_lights[4].directional = false;
		_lights[4].color = vec3(0.0, 0.4, 0.9);
		_lights[4].intensity = 30.0;
*/
    vec3 color = render(ray);
    gl_FragColor = vec4(color, 1.0);

}
