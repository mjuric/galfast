#ifndef __astro__vrml__vrml_h
#define __astro__vrml__vrml_h

#include <astro/util.h>
#include <astro/math/vector.h>

namespace vrml {

std::string header() { return "#VRML V2.0 utf8\n"; }

typedef peyton::math::Vector<float, 3> f3;
typedef peyton::math::Vector<float, 4> f4;
inline OSTREAM(const f3 &v) { return out << v[0] << " " << v[1] << " " << v[2]; }
inline OSTREAM(const f4 &v) { return out << v[0] << " " << v[1] << " " << v[2] << " " << v[3]; }

class Grouping;
class Node {
public:
	Grouping *parent;
public:
	virtual std::ostream &write(std::ostream &out) const = 0;
	virtual std::string type() const = 0;
	virtual ~Node();
	Node(Grouping *parent = NULL);
public: // internal methods
	void setParent(Grouping *parent = NULL);
};
inline OSTREAM(const Node &node) { return node.write(out); }

class Grouping : public Node {
public:
	std::list<Node *> children;
	Grouping(Grouping *parent = NULL) : Node(parent) {}
	virtual ~Grouping();

	Grouping *addChild(Node *node);
	bool removeChild(Node *node);
};

Grouping::~Grouping()
{
	while(children.size())
	{
//		cout << "children [" << children.size() << "] .. deleting " << children.front() << " [" << children.front()->type() << "] \n";
		delete children.front();
	}
}

Grouping *Grouping::addChild(Node *node)
{
	node->setParent(this);
	children.push_back(node);
	return this;
}
bool Grouping::removeChild(Node *node)
{
	std::list<Node *>::iterator child = find(children.begin(), children.end(), node);
	if(child == children.end()) return false;

	children.erase(child);
	// the setting is done directly rather than through
	// setParent because we're either called by setParent (in which case, setParent
	// will set the new parent once we return), or we're called directly, in which
	// case the new parent is NULL
	node->parent = NULL;
	return true;
}

Node::Node(Grouping *p)
{
	parent = NULL;
	if(p != NULL) { p->addChild(this); }
}

// this method MUST NOT be called directly
// use addChild()/removeChild() on the parent instead
// this method is ONLY called from Grouping::addChild
void Node::setParent(Grouping *p)
{
	if(parent != NULL) { parent->removeChild(this); }
	parent = p;
}

Node::~Node()
{
	setParent(NULL);
}

class Transform : public Grouping {
public:
	f3 translation;
	f4 rotation;
	f3 scale;
public:
	virtual std::ostream &write(std::ostream &out) const;
	Transform(Grouping *parent = NULL) : Grouping(parent), translation(0.), rotation(0.), scale(1.) {}

	Transform *setTranslation(float x, float y, float z) { translation = f3(x, y, z); return this; }
	Transform *setScale(float x, float y, float z) { scale = f3(x, y, z); return this; }
	Transform *setRotation(float x, float y, float z, float w) { rotation = f4(x, y, z, w); return this; }
	Transform *setRotation(const f4 &r) { rotation = r; return this; }

	virtual std::string type() const { return "Transform"; }
};
inline Transform *transform() { return (new Transform); }

std::ostream &Transform::write(std::ostream &out) const
{
//	cout << "Transform:\n";
	out << "Transform {\n"
		<< "\ttranslation " << translation << "\n"
		<< "\trotation " << rotation << "\n"
		<< "\tscale " << scale << "\n"
		<< "\tchildren [\n";

//	cout << "\tchildren...\n";
	FOREACH(children) {
//		cout << "\t\tchild... [" << (*i)->type() << "] \n";
		(*i)->write(out);
	}
//	cout << "\tdone...\n";

	out << "\t]"
		<< "}\n";
	return out;
}

class PointLight : public Node {
public:
	f3 location, attenuation, color;
public:
	PointLight(Grouping *parent = NULL)
		: location(0.), attenuation(1., 0, 0), color(1., 1., 1.), Node(parent)
	{}

	PointLight &setLocation(float x, float y, float z) { location = f3(x, y, z); return *this; }

	PointLight *clone(Grouping *parent = NULL) {
		PointLight *pl = new PointLight(parent);
		*pl = *this;
		pl->parent = parent;
		return pl;
	}
	virtual std::ostream &write(std::ostream &out) const;
	virtual std::string type() const { return "PointLight"; }
};

std::ostream &PointLight::write(std::ostream &out) const
{
	out << "PointLight {\n"
		<< "\tattenuation " << attenuation << "\n"
		<< "\tcolor " << color << "\n"
		<< "\tlocation " << location << "\n"
		<< "}\n";
	return out;
}
inline OSTREAM(const PointLight &pl) { return pl.write(out); }

class Shape : public Node
{
public:
	f3 emissiveColor, diffuseColor;
	float transparency;
public:
	Shape(Grouping *parent = NULL)
		: Node(parent), emissiveColor(0, 0, 0), diffuseColor(.8, .8, .8), transparency(0)
	{ }

	Shape *setEmissiveColor(float x, float y, float z) { emissiveColor = f3(x, y, z); return this; }
	Shape *setEmissiveColor(const f3 &ec) { emissiveColor = ec; return this; }
	Shape *setDiffuseColor(float x, float y, float z) { diffuseColor = f3(x, y, z); return this; }
	Shape *setTransparency(float t) { transparency = t; return this; }

	void out_open(std::ostream &out) const
	{
		out << "Shape {\n"
			<< "\tappearance Appearance {\n"
			<< "\t\tmaterial Material {\n"
			<< "\t\t\temissiveColor " << emissiveColor << "\n"
			<< "\t\t\tdiffuseColor " << diffuseColor << "\n"
			<< "\t\t\ttransparency " << transparency << "\n"
			<< "\t\t}\n"
			<< "\t}\n"
			<< "\tgeometry\n";
	}
	void out_close(std::ostream &out) const
	{
		out << "}\n";
	}
};

class Sphere : public Shape
{
public:
	float radius;
public:
	Sphere(Grouping *parent = NULL, double r = 1)
		: Shape(parent), radius(r)
	{}

	virtual std::ostream &write(std::ostream &out) const;
	virtual std::string type() const { return "Sphere"; }
};
inline Sphere *sphere(double r = 1) { return new Sphere(NULL, r); }

std::ostream &Sphere::write(std::ostream &out) const
{
	out_open(out);
	out << "Sphere {\n"
		<< "\tradius " << radius <<"\n"
		<< "}\n";
	out_close(out);
	return out;
}

class Cylinder : public Shape
{
public:
	float radius, height;
public:
	Cylinder(Grouping *parent = NULL, double r = 1, double h = 2)
		: Shape(parent), radius(r), height(h)
	{}

	virtual std::ostream &write(std::ostream &out) const;
	virtual std::string type() const { return "Cylinder"; }
};
inline Cylinder *cylinder(double r = 1, double h = 2) { return new Cylinder(NULL, r, h); }

std::ostream &Cylinder::write(std::ostream &out) const
{
	out_open(out);
	out << "Cylinder {\n"
		<< "\tradius " << radius <<"\n"
		<< "\theight " << height <<"\n"
		<< "}\n";
	out_close(out);
	return out;
}

class Cone : public Shape
{
public:
	float bottomRadius, height;
public:
	Cone(Grouping *parent = NULL, double bottom = 1, double h = 2)
		: Shape(parent), bottomRadius(bottom), height(h)
	{}

	virtual std::ostream &write(std::ostream &out) const;
	virtual std::string type() const { return "Cone"; }
};
inline Cone *cone(double r = 1, double h = 2) { return new Cone(NULL, r, h); }

std::ostream &Cone::write(std::ostream &out) const
{
	out_open(out);
	out << "Cone {\n"
		<< "\tbottomRadius " << bottomRadius <<"\n"
		<< "\theight " << height <<"\n"
		<< "}\n";
	out_close(out);
	return out;
}

class IndexedFaceSet : public Shape
{
public:
	std::vector<f3> vert;
	std::vector<int> idx;
	double epsilon;
	bool solid;
public:
	IndexedFaceSet(Grouping *parent = NULL, double eps = 0.) : Shape(parent), epsilon(eps), solid(true) {}

	void add(std::vector<peyton::math::V3> &poly);
	
	virtual std::ostream &write(std::ostream &out) const;
	virtual std::string type() const { return "IndexedFaceSet"; }
};

void IndexedFaceSet::add(std::vector<peyton::math::V3> &poly)
{
	FOREACH(poly) {
		int j;
		for(j = 0; j != vert.size(); j++) {
			if(epsilon == 0) {
				if(vert[j] == *i) { break; }
			} else {
				if(abs(vert[j] - *i) < epsilon) { break; }
			}
		}
		if(j == vert.size()) { vert.push_back(*i); }
		idx.push_back(j);
	}
	idx.push_back(-1);
}

std::ostream &IndexedFaceSet::write(std::ostream &out) const
{
	out_open(out);

	out << "IndexedFaceSet {\n";
	out << (solid ? "\tsolid TRUE\n" : "\tsolid FALSE\n");
	out << "\tcoord Coordinate {\n";
	out << "\t\tpoint[\n";

	FOREACH(vert) {
		const f3 &v = *i;
		out << "\t\t\t" << v << "\n";
	}
	
	out << "\t\t]\n";
	out << "\t}\n";
	out << "\tcoordIndex [\n";
	bool first = true;
	FOREACH(idx) {
		if(first) { out << "\t\t"; first = false; }
		out << *i;
		if(*i == -1) { first = true; out << "\n"; }
		else { out << " "; }
	}

	out << "\t]\n";
	out << "}\n";

	out_close(out);

	return out;
}

inline OSTREAM(const IndexedFaceSet &ifs) { return ifs.write(out); }

}

#endif
