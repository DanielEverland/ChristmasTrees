#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <unordered_set>
#include <queue>
#include <memory>
#include <climits>
#include <algorithm>

using namespace std;

struct Vertex;

struct Edge
{
	Edge(int weight, int capacity, shared_ptr<Vertex> source, shared_ptr<Vertex> target) : Weight(weight), Capacity(capacity), Source(std::move(source)), Target(std::move(target)) { }

	int Weight;
	int Capacity;
	shared_ptr<Vertex> Source;
	shared_ptr<Vertex> Target;

	void AddWeight(int weight, bool forwardsDirection)
	{
		Weight = forwardsDirection ? Weight + weight : Weight - weight;
	}

	int GetRemainingCapacity(bool forwardsDirection) const
	{
		return forwardsDirection ? Capacity - Weight : Weight;
	}

	bool CanTraverse(bool forwardsDirection) const
	{
		return (!IsFull() && forwardsDirection) || (!IsEmpty() && !forwardsDirection);
	}

	bool IsFull() const
	{
		return Weight >= Capacity;
	}

	bool IsEmpty() const
	{
		return Weight <= 0;
	}
};

struct Vertex
{
	Vertex(string name) : Name(std::move(name)) { }

	string Name;
	vector<shared_ptr<Edge>> OutgoingEdges;
	vector<shared_ptr<Edge>> IngoingEdges;
};

struct PathEntry
{
	PathEntry(bool forward, shared_ptr<PathEntry> from, shared_ptr<Edge> edge) : Forward(forward), From(std::move(from)), PathEdge(std::move(edge)) { }

	bool Forward;
	shared_ptr<PathEntry> From;
	shared_ptr<Edge> PathEdge;

	shared_ptr<Vertex> GetTarget() const
	{
		return Forward ? PathEdge->Target : PathEdge->Source;
	}

	shared_ptr<Vertex> GetSource() const
	{
		return Forward ? PathEdge->Source : PathEdge->Target;
	}
};

struct Hash
{
	size_t operator()(const PathEntry& path) const
	{
		return std::hash<string>{}(path.From->GetTarget()->Name) + std::hash<string>{}(path.From->GetSource()->Name) + std::hash<bool>{}(path.Forward);
	}
};

bool operator==(const PathEntry& a, const PathEntry& b)
{
	return a.PathEdge->Source->Name == b.PathEdge->Source->Name && a.PathEdge->Target->Name == b.PathEdge->Target->Name && a.Forward == b.Forward;
}

struct Path
{
	vector<Edge> SelectedEdges;
	int Weight;
};

void ReadLineInts();
void ReadRoomSize();
void CreateVertices();
void CreateEdges();
void CreateNewEdge(shared_ptr<Vertex> from, shared_ptr<Vertex> to, int capacity);
void CreateNewEdge(int fromIdx, int toIdx, int capacity);
bool FindPath(shared_ptr<PathEntry>& finalPathEntry);


int n, m;
string StringBuffer;
vector<int> IntBuffer;

// 0				= s
// 1 => m			= columns
// m + 1 => m + n	= rows
// m + n			= t
vector<shared_ptr<Vertex>> Vertices;

int main()
{
    ReadRoomSize();
	CreateVertices();
	CreateEdges();

	shared_ptr<PathEntry> path;
	while(FindPath(path))
	{
		shared_ptr<PathEntry> cachedPath = path;
		int minWeight = INT_MAX;
		while(path != nullptr)
		{
			const int remainingCapacity = path->PathEdge->GetRemainingCapacity(path->Forward);
			//cout << "Checking  " << path->PathEdge->Source->Name << " -> " << path->PathEdge->Target->Name << " with remaining capacity " << remainingCapacity << endl;
			if(remainingCapacity < minWeight)
			{
				minWeight = remainingCapacity;
				//cout << "min weight is now " << minWeight << endl;
			}				
				
			path = path->From;
		}

		path = cachedPath;

		while(path != nullptr)
		{
			path->PathEdge->AddWeight(minWeight, path->Forward);
			//cout << "Setting weight of edge " << path->PathEdge->Source->Name << " -> " << path->PathEdge->Target->Name << " to " << path->PathEdge->Weight << endl;
			path = path->From;
		}

		path = nullptr;
	}

	int sum = 0;
	for(shared_ptr<Edge>& edge : Vertices[0]->OutgoingEdges)
	{
		sum += edge->Weight;
	}
	cout << sum << endl;
	return 0;
}

bool FindPath(shared_ptr<PathEntry>& finalPathEntry)
{
	queue<shared_ptr<PathEntry>> openPaths;
	unordered_set<PathEntry, Hash> closed;
	for(shared_ptr<Edge>& edge : Vertices[0]->OutgoingEdges)
	{
		if(edge->CanTraverse(true))
		{
			//cout << "Adding edge which points to " << edge->Target->Name << endl;
			openPaths.emplace(make_shared<PathEntry>(true, nullptr, edge));
		}
	}

	// BFS
	while(!openPaths.empty())
	{
		shared_ptr<PathEntry> currPath = openPaths.front();
		openPaths.pop();

		//cout << "Checking " << currPath->GetSource()->Name << " -> " << currPath->GetTarget()->Name << endl;

		for(shared_ptr<Edge>& edge : currPath->GetTarget()->OutgoingEdges)
		{
			PathEntry path(true, currPath, edge);
			if(closed.find(path) != closed.end())
				continue;

			closed.emplace(path);

			if(edge->Target == Vertices[Vertices.size() - 1] && edge->CanTraverse(true))
			{
				finalPathEntry = make_shared<PathEntry>(true, currPath, edge);
				return true;
			}
			else if(edge->CanTraverse(true))
			{
				openPaths.emplace(make_shared<PathEntry>(path));
			}
		}

		for(shared_ptr<Edge>& edge : currPath->GetTarget()->IngoingEdges)
		{
			PathEntry path(false, currPath, edge);
			if(closed.find(path) != closed.end())
				continue;

			closed.emplace(path);

			if(edge->CanTraverse(false) && edge != currPath->PathEdge)
			{
				openPaths.emplace(make_shared<PathEntry>(path));	
			}
		}
	}

	finalPathEntry = nullptr;
	return false;
}

void CreateEdges()
{
	unordered_set<int> bannedCells;
	for(int i = 0; i < n; i++)
	{
		ReadLineInts();
		bannedCells.clear();
		int max = IntBuffer[0];
		for(int j = 0; j < max; j++)
		{
		 	bannedCells.emplace(IntBuffer[j + 1]);
		}

		for(int k = 0; k < m; k++)
		{
			if(bannedCells.find(k) == bannedCells.end())
			{
				CreateNewEdge(1 + k, m + i + 1, INT_MAX);
			}
		}
	}
}

void CreateVertices()
{
	// Push s to vector
	Vertices.push_back(make_shared<Vertex>("s"));

	int i = 0;
	// All columns
	for(i = 0; i < m; i++)
	{
		Vertices.push_back(make_shared<Vertex>("Column " + to_string(i)));
		CreateNewEdge(0, i + 1, 1);
	}

	// All rows
	shared_ptr<Vertex> t = make_shared<Vertex>("t");
	for(i = 0; i < n; i++)
	{
		Vertices.push_back(make_shared<Vertex>("Row " + to_string(i)));
		CreateNewEdge(Vertices[m + i + 1], t, 2);
	}

	Vertices.push_back(t);
}

void ReadRoomSize()
{
	ReadLineInts();
	n = IntBuffer[0];
	m = IntBuffer[1];
}

void ReadLineInts()
{
	getline(cin, StringBuffer);
	IntBuffer.clear();
	string currentNumberAsString;
	for (const char currentCharacter : StringBuffer)
	{
		if(isspace(currentCharacter) && currentNumberAsString.length() > 0)
		{
			IntBuffer.push_back(atoi(currentNumberAsString.c_str()));
            currentNumberAsString.clear();
		}
        else if(isdigit(currentCharacter))
        {
	        currentNumberAsString += currentCharacter;
        }
	}

	if(currentNumberAsString.length() > 0)
		IntBuffer.push_back(atoi(currentNumberAsString.c_str()));
}

void CreateNewEdge(int fromIdx, int toIdx, int capacity)
{
	//cout << "Creating edge from " << Vertices[fromIdx]->Name << " to " << Vertices[toIdx]->Name << endl;
	CreateNewEdge(Vertices[fromIdx], Vertices[toIdx], capacity);
}

void CreateNewEdge(shared_ptr<Vertex> from, shared_ptr<Vertex> to, int capacity)
{
	shared_ptr<Edge> newEdge = make_shared<Edge>(0, capacity, from, to);

	from->OutgoingEdges.push_back(newEdge);
	to->IngoingEdges.push_back(newEdge);

	//cout << "Created edge from " << newEdge->Source->Name << " to " << newEdge->Target->Name << endl;
}