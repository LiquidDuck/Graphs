#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <algorithm>
#include <stack>
#include <queue>
#include <tuple>

using namespace std;

class DSU
{
private:
	vector<int> p;
	vector<int> rank;

	void makeSet(int oldWeight)
	{
		p[oldWeight] = oldWeight;
		rank[oldWeight] = 0;
	}

public:
	DSU(int N) : p(N + 1), rank(N + 1)
	{
		for (int i = 0; i <= N; ++i) {
			makeSet(i);
		}
	}
	int find(int oldWeight)
	{
		return (oldWeight == p[oldWeight] ? oldWeight : p[oldWeight] = find(p[oldWeight]));
	}
	void unite(int oldWeight, int y)	{		if ((oldWeight = find(oldWeight)) == (y = find(y)))
			return;

		if (rank[oldWeight] <  rank[y])
			p[oldWeight] = y;
		else
			p[y] = oldWeight;

		if (rank[oldWeight] == rank[y])
			++rank[oldWeight];	}	bool checkConnectivity()
	{
		for (int i = 1; i < p.size() - 1; i++)
		{
			if (find(p[i]) != find(p[i + 1]))
			{
				return false;
			}
		}
		return true;
	}
};

class Graph
{
private:
	char type;
	int N = 0;
	int M = 0;
	bool W;
	bool D; 
	vector<vector<int>>  adjMatrix;
	vector<map<int,int>> adjList;
	vector<tuple<int,int,int>> listOfEdges;

	static bool cmpTuple(const tuple<int,int,int> &t1, const tuple<int,int,int> &t2)
	{
		if (get<2>(t1) != get<2>(t2))
		{
			return (get<2>(t1) < get<2>(t2));
		}
		else
			if (get<0>(t1) != get<0>(t2))
			{
				return (get<0>(t1) < get<0>(t2));
			}
			else
			{
				return (get<1>(t1) < get<1>(t2));
			}
	}
	bool checkBridge(int n1, int n2)
	{
		Graph graph('L', D, W, N, M);
		DSU dsu(N);
		transformToAdjList();

		for (int i = 0; i < N; i++)
		{
			graph.adjList[i] = adjList[i];
		}
		graph.removeEdge(n1, n2);

		for (int i = 0; i < N; i++)
		{
			int a = i + 1;
			for (auto j = graph.adjList[i].begin(); j != graph.adjList[i].end(); j++)
			{
				int b = j->first;
				if (dsu.find(a) != dsu.find(b))
				{
					dsu.unite(a, b);
				}
			}
		}
		return !dsu.checkConnectivity();
	}
	bool kuhnDFS(int v, vector<bool>& used, vector<int>& curBipart) 
	{
		if (used[v])
		{
			return false;
		}
		used[v] = true;

		for (auto& i : adjList[v - 1]) 
		{
			int to = i.first;
			if (curBipart[to] == -1 || kuhnDFS(curBipart[to], used, curBipart)) 
			{
				curBipart[to] = v;
				return true;
			}
		}
		return false;
	}

public: 
	Graph() {}
	Graph(char ntype, bool nD, bool nW, int nN, int nM) : type(ntype), D(nD), W(nW), N(nN), M(nM)
	{
		switch (type)
		{
		case 'C':
		{
			adjMatrix.resize(N);
			for (int i = 0; i < N; i++)
				adjMatrix[i].resize(N);
			break;
		}
		case 'L':
		{
			adjList.resize(N);
			break;
		}
		case 'E':
		{
			listOfEdges.resize(M);
			break;
		}
		}
	}
	Graph(int nN):type('E'), D(0), W(1), N(nN), M(0) {}
	Graph(const Graph& graph) : type(graph.type), D(graph.D), W(graph.W), N(graph.N), M(graph.M), adjMatrix(graph.adjMatrix),
		adjList(graph.adjList), listOfEdges(graph.listOfEdges) {}
	~Graph() {}

	void readGraph(string fileName)
	{
		adjList.clear();
		adjMatrix.clear();
		listOfEdges.clear();
		ifstream inS(fileName);
		inS >> type;
		
		switch (type)
		{
		case 'C':
		{
			inS >> N >> D >> W;
			adjMatrix.resize(N);
			for (int i = 0; i < N; i++)
				adjMatrix[i].resize(N);

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					int w;
					inS >> w;
					w ? M++ : 1;
					adjMatrix[i][j] = w;
				}
			}
			if (!D)
				M /= 2;
			break;
		}
		case 'L':
		{
			inS >> N >> D >> W;
			adjList.resize(N);
			string s;
			getline(inS, s);

			for (int i = 0; i < N; i++)
			{
				getline(inS, s);
				istringstream iss(s);
				int b, w = 1;
				if (s.size() != 0)
				{
					while (iss && s.size() != 0)
					{
						iss >> b;
						if (W)
						{
							iss >> w;
						}
						adjList[i][b] = w;
					}
				}
				M += adjList[i].size();
			}
			break;
		}
		case 'E':
		{
			inS >> N >> M >> D >> W;
			listOfEdges.resize(M);

			for (int i = 0; i < M; i++)
			{
				int a, b, w = 1;
				inS >> a >> b;
				if (W)
				{
					inS >> w;
				}
				get<0>(listOfEdges[i]) = a;
				get<1>(listOfEdges[i]) = b;
				get<2>(listOfEdges[i]) = w;
			}
			break;
		}
		}

		inS.close();
	}
	void addEdge(int from, int to, int weight)
	{
		if (!D && from > to)
		{
			swap(from, to);
		}

		switch (type)
		{
		case 'C':
		{
			M++;
			adjMatrix[from - 1][to - 1] = weight;
			D ? 1 : adjMatrix[to - 1][from - 1] = weight;
			break;
		}
		case 'L':
		{
			M++;
			adjList[from - 1][to] = weight;
			D ? 1 : adjList[to - 1][from] = weight;
			break;
		}
		case 'E':
		{
			M++;
			listOfEdges.resize(M);
			get<0>(listOfEdges[M - 1]) = from;
			get<1>(listOfEdges[M - 1]) = to;
			get<2>(listOfEdges[M - 1]) = weight;
			break;
		}
		}
	}
	void removeEdge(int from, int to)
	{
		switch (type)
		{
		case 'C':
		{
			if (adjMatrix[from - 1][to - 1] != 0)
			{
				adjMatrix[from - 1][to - 1] = 0;
				D ? 1 : adjMatrix[to - 1][from - 1] = 0;
				M--;
			}
			break;
		}
		case 'L':
		{
			adjList[from - 1].erase(to);
			M--;
			D ? 1 : adjList[to - 1].erase(from);
			break;
		}
		case 'E':
		{
			if (!D && from > to)
			{
				swap(from, to);
			}

			auto edge = find_if(listOfEdges.begin(), listOfEdges.end(), [from, to](tuple<int,int,int>& edge) 
			{
				return (std::get<0>(edge) == from) && (std::get<1>(edge) == to);
			});

			listOfEdges.erase(edge);
			M--;
			break;
		}
		}
	}
	int changeEdge(int from, int to, int newWeight)
	{
		int oldWeight = 0;
		switch (type)
		{
		case 'C':
		{
			oldWeight = adjMatrix[from - 1][to - 1];
			adjMatrix[from - 1][to - 1] = newWeight;
			D ? 1 : adjMatrix[to - 1][from - 1] = newWeight;
			break;
		}
		case 'L':
		{
			oldWeight = adjList[from - 1][to];
			adjList[from - 1][to] = newWeight;
			D ? 1 : adjList[to - 1][from] = newWeight;
			break;
		}
		case 'E':
		{
			if (!D && from > to)
			{
				swap(from, to);
			}

			auto edge = find_if(listOfEdges.begin(), listOfEdges.end(), [from, to](tuple<int,int,int>& edge)
			{
				return (get<0>(edge) == from) && (get<1>(edge) == to);
			});

			oldWeight = get<2>(*edge);
			get<2>(*edge) = newWeight;
			break;
		}
		}
		return oldWeight;
	}	
	void transformToAdjMatrix()
	{
		switch (type)
		{
		case 'C':
			break;
		case 'L':
		{
			adjMatrix.resize(N);
			for (int i = 0; i < N; i++)
			{
				adjMatrix[i].resize(N);
			}

			for (int i = 0; i < N; i++)
			{
				for (auto it = adjList[i].begin(); it != adjList[i].end(); it++)
				{
					int b = it->first;
					int w = it->second;
					adjMatrix[i][b - 1] = w;
					D ? 1 : adjMatrix[b - 1][i] = w;
				}
			}

			adjList.clear();
			adjList.shrink_to_fit();
			type = 'C';
			break;
		}
		case 'E':
		{
			adjMatrix.resize(N);
			for (int i = 0; i < N; i++)
			{
				adjMatrix[i].resize(N);
			}

			for (int i = 0; i < M; i++)
			{
				int a = get<0>(listOfEdges[i]);
				int b = get<1>(listOfEdges[i]);
				int w = get<2>(listOfEdges[i]);

				adjMatrix[a - 1][b - 1] = w;
				D ? 1 : adjMatrix[b - 1][a - 1] = w;
			}

			listOfEdges.clear();
			listOfEdges.shrink_to_fit();
			type = 'C';
			break;
		}
		}
	}
	void transformToAdjList()
	{
		switch (type)
		{
		case 'C':
		{
			adjList.resize(N);
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (adjMatrix[i][j] != 0)
					{
						adjList[i][j + 1] = adjMatrix[i][j];
						D ? 1 : adjList[j][i + 1] = adjMatrix[i][j];
					}
				}

			adjMatrix.clear();
			adjMatrix.shrink_to_fit();
			type = 'L';
			break;
		}
		case 'L':
			break;
		case 'E':
		{
			adjList.resize(N);

			for (int i = 0; i < M; i++)
			{
				int a = get<0>(listOfEdges[i]);
				int b = get<1>(listOfEdges[i]);
				int w = get<2>(listOfEdges[i]);

				adjList[a - 1][b] = w;
				D ? 1 : adjList[b - 1][a] = w;
			}

			listOfEdges.clear();
			listOfEdges.shrink_to_fit();
			type = 'L';
			break;
		}
		}
	}
	void transformToListOfEdges()
	{
		switch (type)
		{
		case 'C':
		{
			listOfEdges.resize(M);
			int k = 0;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (adjMatrix[i][j] != 0)
					{
						if (!D && i < j || D)
						{
							get<0>(listOfEdges[k]) = i + 1;
							get<1>(listOfEdges[k]) = j + 1;
							get<2>(listOfEdges[k]) = adjMatrix[i][j];
							k++;
						}
					}
				}

			adjMatrix.clear();
			adjMatrix.shrink_to_fit();
			type = 'E';
			break;
		}
		case 'L':
		{
			listOfEdges.resize(M);
			int k = 0;

			for (int i = 0; i < N; i++)
			{
				for (auto it = adjList[i].begin(); it != adjList[i].end(); it++)
				{
					int b = it->first;
					int w = it->second;
					if (!D && i < b || D)
					{
						get<0>(listOfEdges[k]) = i + 1;
						get<1>(listOfEdges[k]) = b;
						get<2>(listOfEdges[k]) = w;
						k++;
					}
				}
			}

			adjList.clear();
			adjList.shrink_to_fit();
			type = 'E';
			break;
		}
		case 'E':
			break;
		}
	}
	void writeGraph(string fileName)
	{
		ofstream outS(fileName);

		switch (type)
		{
		case 'C':
		{
			outS << type << ' ' << N << endl << D << ' ' << W << endl;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N - 1; j++)
				{
					outS << adjMatrix[i][j] << ' ';
				}
				outS << adjMatrix[i][N - 1] << endl;
			}
			break;
		}
		case 'L':
		{
			outS << type << ' ' << N << endl << D << ' ' << W << endl;

			for (int i = 0; i < N; i++)
			{
				int adjListSize = adjList[i].size();
				int k = 0;

				for (auto it = adjList[i].begin(); it != adjList[i].end(); it++)
				{
					if (k == adjListSize - 1)
					{
						outS << it->first;
						if (W)
						{
							outS << ' ' << it->second;
						}
					}
					else
					{
						outS << it->first << ' ';
						if (W)
							outS << it->second << ' ';
					}
					k++;
				}
				outS << endl;
			}
			break;
		}
		case 'E':
		{
			outS << type << ' ' << N << ' ' << M << endl << D << ' ' << W << endl;

			for (int i = 0; i < M; i++)
			{
				outS << get<0>(listOfEdges[i]) << ' ' << get<1>(listOfEdges[i]);
				if (W)
				{
					outS << ' ' << get<2>(listOfEdges[i]);
				}
				outS << endl;
			}
			break;
		}
		}

		outS.close();
	}

	Graph getSpanningTreePrima()
	{
		Graph graph(N);
		vector<int> key(N + 1, 1000000);
		vector<int> parent(N + 1, -1);
		set<int> used;
		set<pair<int,int>> queue;

		queue.insert(make_pair(0, 1));
		key[1] = 0;
		transformToAdjList();

		for (int i = 1; i <= N; i++)
		{
			used.insert(i);
		}

		for (int i = 0; i < N; i++)
		{
			if (queue.empty())
			{
				queue.insert(make_pair(0, *used.begin()));
			}
			int v = queue.begin()->second;
			queue.erase(queue.begin());

			for (auto j = adjList[v - 1].begin(); j != adjList[v - 1].end(); j++)
			{
				int to = j->first;
				int cost = j->second;
				if (cost < key[to] && used.find(to) != used.end())
				{
					queue.erase(make_pair(key[to], to));
					key[to] = cost;
					parent[to] = v;
					queue.insert(make_pair(key[to], to));
				}
			}

			used.erase(v);
			if (!queue.empty())
			{
				graph.addEdge(parent[queue.begin()->second], queue.begin()->second, queue.begin()->first);
			}
		}
		return graph;
	}
	Graph getSpanningTreeKruscal()
	{
		Graph graph(N);
		transformToListOfEdges();
		sort(listOfEdges.begin(), listOfEdges.end(), cmpTuple);
		DSU dsu(N);

		for (int i = 0; i < M; i++)
		{
			int a = get<0>(listOfEdges[i]);
			int b = get<1>(listOfEdges[i]);
			int w = get<2>(listOfEdges[i]);
			if (dsu.find(a) != dsu.find(b))
			{
				graph.addEdge(a, b, w);
				dsu.unite(a, b);
			}
		}
		return graph;
	}
	Graph getSpanningTreeBoruvka()	{		Graph graph(N);
		transformToListOfEdges();
		DSU dsu(N);
		int k = N;

		while (k > 1)
		{
			bool flag = true;
			vector<int> key(N + 1, -1);

			for (int i = 0; i < M; i++)
			{
				int a = dsu.find(get<0>(listOfEdges[i]));
				int b = dsu.find(get<1>(listOfEdges[i]));
				int w = get<2>(listOfEdges[i]);

				if (a != b)
				{
					if (key[a] == -1 || w < get<2>(listOfEdges[key[a]]))
					{
						key[a] = i;
						flag = false;
					}
					if (key[b] == -1 || w < get<2>(listOfEdges[key[b]]))
					{
						key[b] = i;
						flag = false;
					}
				}
			}

			if (flag) 
				break;

			for (int i = 1; i <= N; i++)
			{
				if (key[i] != -1)
				{
					int a = get<0>(listOfEdges[key[i]]);
					int b = get<1>(listOfEdges[key[i]]);
					int w = get<2>(listOfEdges[key[i]]);
					if (dsu.find(a) != dsu.find(b))
					{
						graph.addEdge(a, b, w);
						dsu.unite(a, b);
						k--;
					}
				}
			}
		}
		return graph;	}

	int checkEuler(bool &circleExist)
	{
		DSU dsu(N);
		bool dsuFlag = false;
		int res = 1, odd = 0;
		vector<int> power(N + 1, 0);
		transformToListOfEdges();

		for (int i = 0; i < M; i++)
		{
			power[get<0>(listOfEdges[i])]++;
			power[get<1>(listOfEdges[i])]++;
		}

		for (int i = 1; i <= N; i++)
		{
			if (power[i] % 2 == 1)
			{
				odd++;
				res = i;
			}
		}

		for (int i = 0; i < M; i++)
		{
			int a = get<0>(listOfEdges[i]);
			int b = get<1>(listOfEdges[i]);
			int w = get<2>(listOfEdges[i]);
			if (dsu.find(a) != dsu.find(b))
			{
				dsu.unite(a, b);
			}
		}

		dsuFlag = dsu.checkConnectivity();
		odd == 0 ? circleExist = true : circleExist = false;
		if ((odd == 0 || odd == 2) && dsuFlag)
		{
			return res;
		}
		else
		{
			return 0;
		}
	}
	vector<int> getEuleranTourFleri()
	{
		Graph graph('E', D, W, N, M);
		vector<int> res(M + 1);
		vector<int> power(N + 1, 0);
		transformToListOfEdges();

		for (int i = 0; i < M; i++)
		{
			power[get<0>(listOfEdges[i])]++;
			power[get<1>(listOfEdges[i])]++;
			get<0>(graph.listOfEdges[i]) = get<0>(listOfEdges[i]);
			get<1>(graph.listOfEdges[i]) = get<1>(listOfEdges[i]);
			get<2>(graph.listOfEdges[i]) = get<2>(listOfEdges[i]);
		}

		bool circleFlag;
		int start = checkEuler(circleFlag);
		res[0] = start;
		transformToAdjList();
		graph.transformToAdjList();
		int current = start;

		for (int i = 0; i < M; i++)
		{
			bool stepFlag = false;
			for (auto j = graph.adjList[current - 1].begin(); j != graph.adjList[current - 1].end(); j++)
			{
				int neigh = j->first;
				if (!graph.checkBridge(current, neigh) || power[current] == 1)
				{
					stepFlag = true;
					power[current]--;
					graph.removeEdge(current, neigh);
					current = neigh;
					res[i + 1] = current;
					power[neigh]--;
					break;
				}
			}

			if (!stepFlag) 
				for (auto j = graph.adjList[current - 1].begin(); j != graph.adjList[current - 1].end(); j++)
				{
					int neigh = j->first;
					if (!graph.checkBridge(current, neigh) || power[neigh] != 1)
					{
						power[current]--;
						graph.removeEdge(current, neigh);
						current = neigh;
						res[i + 1] = current;
						power[neigh]--;
						break;
					}
				}
			stepFlag = false;
		}
		return res;
	}
	vector<int> getEuleranTourEffective()	{		vector<int> result;
		bool circleExist;
		int start = checkEuler(circleExist);
		if (start == 0) 
		{
			return result;
		}

		transformToAdjList();
		Graph graph(*this);
		stack<int> stack;
		stack.push(start);

		while (!stack.empty())
		{
			int curV = stack.top();
			for (auto& edge : adjList[curV - 1]) 
			{
				int nextV = edge.first;
				stack.push(nextV);
				graph.removeEdge(curV, nextV);
				break;
			}
			if (curV == stack.top()) 
			{
				stack.pop();
				result.push_back(curV);
			}
		}
		return result;	}

	int checkBipart(vector<char> &marks)
	{
		transformToAdjList();
		queue<int> queue;

		if (N < 2)
		{
			return false;
		}

		for (int i = 1; i <= N; i++)
		{
			if (marks[i] != 'A' && marks[i] != 'B') 
			{
				marks[i] = 'A';
				queue.push(i);
				while (!queue.empty()) 
				{
					int prevV = queue.front();
					queue.pop();

					for (auto& i : adjList[prevV - 1]) 
					{
						if (marks[i.first] != 'A' && marks[i.first] != 'B') 
						{
							marks[i.first] = marks[prevV] == 'A' ? 'B' : 'A';
							queue.push(i.first);
						}
						else 
							if (marks[i.first] == marks[prevV])
							{
								return false;
							}
					}
				}
			}
		}
		return true;
	}
	vector<pair<int,int>> getMaximumMatchingBipart()	{
		transformToAdjList();
		vector<pair<int,int>> res;
		vector<int> curB(N + 1, -1);
		vector<char> marks(N + 1, 'Z');
		vector<bool> used(N + 1, false);
		vector<bool> used1(N + 1);

		for (int i = 1; i <= N; i++) 
		{
			if (marks[i] == 'B')
			{
				continue;
			}
			for (auto& j : adjList[i - 1]) 
			{
				if (curB[j.first] == -1)
				{
					curB[j.first] = i;
					used1[i] = true;
					break;
				}
			}
		}

		if (checkBipart(marks)) 
		{
			for (int i = 1; i <= N; i++) 
			{
				if (used1[i] || marks[i] == 'B')
				{
					continue;
				}
				used.assign(N + 1, false);
				kuhnDFS(i, used, curB);
			}
			for (int i = 1; i <= N; i++)
			{
				if (marks[i] == 'B' && curB[i] != -1)
				{
					res.emplace_back(curB[i], i);
				}
			}
		}
		return res;
	}
};
