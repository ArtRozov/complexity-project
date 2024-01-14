#include <algorithm>
#include <iostream>
#include <vector>
#include <optional>
#include <map>
#include <stack>

class DSU {
 public:
  DSU() = default;

  DSU(long long vertices_num) {
    dsu_.resize(vertices_num);
    for (long long i = 0; i < vertices_num; ++i) {
      dsu_[i] = i;
    }
  }

  long long GetParent(long long vert) {
    return (vert == dsu_[vert]) ? vert : dsu_[vert] = GetParent(dsu_[vert]);
  }

  bool EqualParent(long long vert1, long long vert2) { return GetParent(vert1) == GetParent(vert2); }


  void Unite(long long vert1, long long vert2) {
    vert1 = GetParent(vert1);
    vert2 = GetParent(vert2);
    dsu_[vert1] = vert2;
  }

  void Reset() {
    for (long long i = 0; i < dsu_.size(); ++i) {
      dsu_[i] = i;
    }
  }

 private:
  std::vector<long long> dsu_;
};

template <bool is_oriented>
class Graph {
 public:
  struct Edge {
    bool operator<(const Edge& other) const { return weight < other.weight; }

    Edge() = default;
    Edge(long long in, long long out) { from = in; to = out; weight = 1; }
    Edge(long long in, long long out, long long w) { from = in; to = out; weight = w; }

    long long from;
    long long to;
    long long weight;
  };

  Graph() = default;

  Graph(size_t vertices_num) {
    adjacency_list_.resize(vertices_num);
    mark_ = std::vector<bool>(vertices_num, false);
  }

  void ResetSettings() {
    max_degree_ = 0;
    for (const auto& edges : adjacency_list_) {
      max_degree_ = std::max<long long>(max_degree_, edges.size());
    }
    mark_ = std::vector<bool>(adjacency_list_.size(), false);
  }

  void ReadGraph(std::istream& in = std::cin) {
    size_t edges_num;
    in >> edges_num;
    for (size_t i = 0; i < edges_num; ++i) {
      long long from;
      long long to;
      long long weight = 1;
      in >> from >> to;
      --from;
      --to;
      AddEdge(from, to, weight);
    }
  }

  void AddEdge(long long from, long long to, long long weight = 1) {
    adjacency_list_[from].push_back(std::pair(weight, to));
    if constexpr (!is_oriented) {
      adjacency_list_[to].push_back(std::pair(weight, from));
    }
  }

  void RemoveEdge(long long from, long long to) {
    adjacency_list_[from].erase(
        std::remove_if(
            adjacency_list_[from].begin(),
            adjacency_list_[from].end(),
            [to](const std::pair<long long, long long>& edge) { return edge.second == to; }
        ),
        adjacency_list_[from].end()
    );

    if constexpr (!is_oriented) {
        adjacency_list_[to].erase(
            std::remove_if(
                adjacency_list_[to].begin(),
                adjacency_list_[to].end(),
                [from](const std::pair<long long, long long>& edge) { return edge.second == from; }
            ),
            adjacency_list_[to].end()
        );
    }
}



  std::vector<Edge> GetEdges() const {
    std::vector<Edge> edges;
    for (size_t i = 0; i < adjacency_list_.size(); ++i) {
      for (size_t j = 0; j < adjacency_list_[i].size(); ++j) {
        edges.push_back(Edge(static_cast<long long>(i),
                         adjacency_list_[i][j].second));
      }
    }
    return edges;
  }

  size_t GetVerticesNum() const { return adjacency_list_.size(); }

  void UpdateMSTWeight(long long cost) { mst_cost_ += cost; }

  long long GetMSTWeight() const { return mst_cost_; }

  size_t GetMaxDegree() const { return max_degree_; }

  size_t GetDegree(long long vert) const { return adjacency_list_[vert].size(); }

  std::vector<size_t> FindPathDFS(size_t start, size_t end) const {
        std::vector<bool> visited(adjacency_list_.size(), false);
        std::vector<size_t> path;
        std::stack<size_t> stack;
        stack.push(start);
        visited[start] = true;
        while (!stack.empty()) {
            size_t current = stack.top();
            stack.pop();
            path.push_back(current);
            if (current == end) { return path; }
            for (auto [w, neighbor] : adjacency_list_[current]) {
                if (!visited[neighbor]) {
                    stack.push(neighbor);
                    visited[neighbor] = true;
                }
            }
        }
        return {};
  }
  
  void MarkComponents() {
    for (long long vert = 0; vert < adjacency_list_.size(); ++vert) {
      mark_[vert] = adjacency_list_[vert].size() + 1 < max_degree_;
    }
  }

  bool IsVertexGood(long long vert) const { return mark_[vert]; }

  void MarkVertexGood(long long vert) { mark_[vert] = true; }

  std::vector<bool> Marks() { return mark_; }

  void ComponentsDivision(DSU& dsu) {
    for (size_t vert1 = 0; vert1 < adjacency_list_.size(); ++vert1) {
      if (!IsVertexGood(vert1)) continue;
      for (auto edge : adjacency_list_[vert1]) {
        auto vert2 = edge.second;
        if (IsVertexGood(vert2) && !dsu.EqualParent(vert1, vert2)) {
          dsu.Unite(vert1, vert2);
        }
      }
    }
  }

  void PrintGraph() {
    size_t vertices_num = adjacency_list_.size();

    std::cout << "Vertices number - " << vertices_num << "\n\n";

    std::cout << "Adjacency Matrix:\n";
    for (size_t i = 0; i < vertices_num; ++i) {
        for (size_t j = 0; j < vertices_num; ++j) {
            auto it = std::find_if(adjacency_list_[i].begin(), adjacency_list_[i].end(),
                                   [j](const auto& edge) { return edge.second == static_cast<long long>(j); });

            if (it != adjacency_list_[i].end()) {
                std::cout << it->first << "\t";
            } else {
                std::cout << "0\t";
            }
        }
        std::cout << "\n";
    }
}


private:
  long long mst_cost_ = 0;
  std::vector<std::vector<std::pair<long long, long long>>> adjacency_list_;
  size_t max_degree_;
  std::vector<bool> mark_;
};

class MSTBuilder {
 public:
  MSTBuilder() = default;

  virtual Graph<false> MstBuilder(const Graph<false>& graph) = 0;
};

class MSTFinder : MSTBuilder {
private:
  struct ImprovementObject {
    long long delete_edge_vert_;
    Graph<false>::Edge improvement_edge_;

    ImprovementObject() = default;

    ImprovementObject(long long to, Graph<false>::Edge improvement_edge) {
      delete_edge_vert_ = to;
      improvement_edge_ = improvement_edge;
    }
  };

public:
  MSTFinder(long long vertices_num) { 
    vertices_num_ = vertices_num;
    dsu_ = DSU(vertices_num);
  }

  Graph<false> MstBuilder(const Graph<false>& graph) override {
    Graph<false> return_graph(graph.GetVerticesNum());
    std::vector<Graph<false>::Edge> edges = graph.GetEdges();
    std::sort(edges.begin(), edges.end());
    for (size_t i = 0; i < edges.size(); ++i) {
      long long from = edges[i].from;
      long long to = edges[i].to;
      long long cost = edges[i].weight;
      if (dsu_.GetParent(from) != dsu_.GetParent(to)) {
        return_graph.AddEdge(from, to, cost);
        return_graph.UpdateMSTWeight(cost);
        dsu_.Unite(from, to);
      }
    }
    return return_graph;
  }

  void MstReInit(Graph<false>& mst) {
    dsu_.Reset();
    mst.ResetSettings();

    mst.MarkComponents();
    mst.ComponentsDivision(dsu_);
  }

  void Improve(Graph<false>& mst, long long vert) {
    auto it = improvements.find(vert);
    if (it != improvements.end()) {
        const auto& obj = it->second;
        if (mst.GetDegree(obj.improvement_edge_.from) == static_cast<long long>(mst.GetMaxDegree()) - 1) {
            Improve(mst, obj.improvement_edge_.from);
        }
        if (mst.GetDegree(obj.improvement_edge_.to) == static_cast<long long>(mst.GetMaxDegree()) - 1) {
            Improve(mst, obj.improvement_edge_.to);
        }
        mst.RemoveEdge(vert, obj.delete_edge_vert_);
        mst.AddEdge(obj.improvement_edge_.from, obj.improvement_edge_.to);
        improvements.erase(it);
    }
}


  bool TryImprove(Graph<false>& mst, long long vert1, long long vert2) {
    auto cycle = mst.FindPathDFS(vert1, vert2);
    for (size_t i = 0; i < cycle.size(); ++i) {
      size_t vert = cycle[i];
      if (!mst.IsVertexGood(vert) && mst.GetDegree(vert) == mst.GetMaxDegree()) {
        ImprovementObject obj(cycle[i + 1 < cycle.size() ? i + 1 : i - 1], Graph<false>::Edge(vert1, vert2));
        improvements[vert] = obj;
        Improve(mst, vert);
        return true;
      }
    }
    // max_degree - 1    process
    if (!cycle.empty()) cycle.pop_back();
    for (size_t i = 0; i < cycle.size(); ++i) {
      size_t vert = cycle[i];
      if (!mst.IsVertexGood(vert) && mst.GetDegree(vert) == static_cast<long long>(mst.GetMaxDegree()) - 1) {
        ImprovementObject obj(cycle[i + 1 < cycle.size() ? i + 1 : i - 1], Graph<false>::Edge(vert1, vert2));
        improvements[vert] = obj;
      }
      mst.MarkVertexGood(vert);
      if (!dsu_.EqualParent(vert1, vert2)) {
        dsu_.Unite(vert, vert2);
      }
    }
    return false;
  }

  bool Equal(long long vert1, long long vert2) { return dsu_.EqualParent(vert1, vert2); }

 private:
  size_t vertices_num_;
  DSU dsu_;
  std::map<size_t, ImprovementObject> improvements;
};

std::pair<std::pair<Graph<false>, Graph<false>>, MSTFinder> MstInit() {
  long long vertices_num;
  std::cin >> vertices_num;
  Graph<false> graph(vertices_num);
  graph.ReadGraph();
  MSTFinder kruskal(vertices_num);
  Graph<false> mst = kruskal.MstBuilder(graph);
  return std::make_pair(std::make_pair(graph, mst), kruskal);
}

void PhasesInit(const Graph<false>& graph, Graph<false>& mst, MSTFinder& mst_finder) {
  bool no_improvement_indicated = false;
  while (!no_improvement_indicated) {  
    mst_finder.MstReInit(mst);
    bool improved = false;
    auto edges = graph.GetEdges();
    bool improvement_edge_found = false;
    while (!improved) {
      bool improvement_edge_found = false;
      for (auto edge : edges) {
        long long vert1 = edge.from;
        long long vert2 = edge.to;
        if (mst.IsVertexGood(vert1) && mst.IsVertexGood(vert2) && !mst_finder.Equal(vert1, vert2)) {
          improvement_edge_found = true;
          if ((improved = mst_finder.TryImprove(mst, vert1, vert2))) {
            break;
          }
        }
      }
      if (!improvement_edge_found) { no_improvement_indicated = true; break; }
    } 
  }
  std::cout << "Maximum degree is " << mst.GetMaxDegree() << '\n';
  mst.PrintGraph();
}

int main() {
  auto [graphes, kruskal] = MstInit();
  auto [graph, mst] = graphes;
  PhasesInit(graph, mst, kruskal);
  return 0;
}
