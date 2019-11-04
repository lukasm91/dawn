#pragma once

#include <cmath>
#include <iostream>
#include <vector>

namespace wstd {
struct identity {
  template <class T>
  constexpr T&& operator()(T&& t) const noexcept {
    return std::forward<T>(t);
  }
};
} // namespace wstd

namespace mylib {
class Vertex;
class Edge;
class Face;

//   .---.---.---.
//   |\ 1|\ 3|\ 5|
//   | \ | \ | \ |
//   |0 \|2 \|4 \|
//   ----'---'---'
enum face_color { upward = 0, downward = 1 };
enum edge_color { horizontal = 0, diagonal = 1, vertical = 2 };

class Vertex {
public:
  Vertex() = default;
  Vertex(int x, int y, int z, int id) : x_(x), y_(y), z_(z), id_(id) {}

  int x() const { return x_; }
  int y() const { return y_; }
  int z() const { return z_; }
  int id() const { return id_; }

  Edge const& edge(size_t i) const;
  Face const& face(size_t i) const;
  Vertex const& vertex(size_t i) const;
  auto edges() const { return edges_; }
  auto faces() const { return faces_; }
  std::vector<const Vertex*> vertices() const;

  void add_edge(Edge& e);
  void add_face(Face& f) { faces_.push_back(&f); }

private:
  int id_;
  int x_;
  int y_;
  int z_;

  std::vector<Edge*> edges_;
  std::vector<Face*> faces_;
};
class Face {
public:
  Face() = default;
  Face(int id, face_color color) : id_(id), color_(color) {}

  int id() const { return id_; }
  face_color color() const { return color_; }

  Vertex const& vertex(size_t i) const;
  Edge const& edge(size_t i) const;
  Face const& face(size_t i) const;
  auto vertices() const { return vertices_; }
  auto edges() const { return edges_; }
  std::vector<const Face*> faces() const;

  void add_edge(Edge& e) { edges_.push_back(&e); }
  void add_vertex(Vertex& v) { vertices_.push_back(&v); }

private:
  int id_;
  face_color color_;

  std::vector<Edge*> edges_;
  std::vector<Vertex*> vertices_;
};
class Edge {
public:
  Edge() = default;
  Edge(int id, edge_color color) : id_(id), color_(color) {}

  int id() const { return id_; }
  edge_color color() const { return color_; }

  Vertex const& vertex(size_t i) const;
  Face const& face(size_t i) const;
  auto faces() const { return faces_; }
  auto vertices() const { return vertices_; }

  void add_vertex(Vertex& v) { vertices_.push_back(&v); }
  void add_face(Face& f) { faces_.push_back(&f); }

  operator bool() const { return id_ >= 0; }

private:
  int id_ = -1;
  edge_color color_;

  std::vector<Vertex*> vertices_;
  std::vector<Face*> faces_;
};

class Grid {
public:
  Grid(int nx, int ny, int nz, bool periodic = false)
      : faces_(2 * nx * ny * nz), vertices_(periodic ? nx * ny * nz : (nx + 1) * (ny + 1) * nz),
        edges_(periodic ? 3 * nx * ny * nz : 3 * (nx + 1) * (ny + 1) * nz), nx_(nx), ny_(ny),
        nz_(nz) {
    auto edge_at = [&](int i, int j, int k, int c) -> Edge& {
      if(periodic)
        return edges_.at(((((j + ny) % ny) * nx + ((i + nx) % nx)) * 3 + c) * nz + k);
      else
        return edges_.at(((j * (nx + 1) + i) * 3 + c) * nz + k);
    };
    auto vertex_at = [&](int i, int j, int k) -> Vertex& {
      if(periodic)
        return vertices_.at((((j + ny) % ny) * nx + ((i + nx) % nx)) * nz + k);
      else
        return vertices_.at((j * (nx + 1) + i) * nz + k);
    };
    auto face_at = [&](int i, int j, int k, int c) -> Face& {
      if(periodic)
        return faces_.at(((((j + ny) % ny) * nx + ((i + nx) % nx)) * 2 + c) * nz + k);
      else
        return faces_.at(((j * nx + i) * 2 + c) * nz + k);
    };

    for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
        for(int i = 0; i < nx; ++i)
          for(int c = 0; c < 2; ++c) {
            auto& f = face_at(i, j, k, c);
            f = Face(&f - faces_.data(), (face_color)c);
          }
    for(int k = 0; k < nz; ++k)
      for(int j = 0; j < (periodic ? ny : ny + 1); ++j)
        for(int i = 0; i < (periodic ? nx : nx + 1); ++i) {
          auto& v = vertex_at(i, j, k);
          v = Vertex(i, j, k, &v - vertices_.data());
        }

    for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
        for(int i = 0; i < nx; ++i) {
          auto& f_uw = face_at(i, j, k, face_color::upward);
          //   .
          //   |\
        //  0| \ 1
          //   |  \ 
        //   ----'
          //     2
          f_uw.add_edge(edge_at(i, j, k, edge_color::vertical));
          f_uw.add_edge(edge_at(i, j, k, edge_color::diagonal));
          f_uw.add_edge(edge_at(i, j + 1, k, edge_color::horizontal));

          //   0
          //   .
          //   |\ 
        //   | \ 
        //   |  \ 
        //   ----' 1
          //  2
          f_uw.add_vertex(vertex_at(i, j, k));
          f_uw.add_vertex(vertex_at(i + 1, j + 1, k));
          f_uw.add_vertex(vertex_at(i, j + 1, k));

          // downward
          auto& f_dw = face_at(i, j, k, face_color::downward);
          //     1
          //   ----
          //   \  |
          //  0 \ |2
          //     \|
          //      ^
          f_dw.add_edge(edge_at(i, j, k, edge_color::diagonal));
          f_dw.add_edge(edge_at(i, j, k, edge_color::horizontal));
          f_dw.add_edge(edge_at(i + 1, j, k, edge_color::vertical));

          //        1
          // 0 ----
          //   \  |
          //    \ |
          //     \|
          //      ^ 2
          f_dw.add_vertex(vertex_at(i, j, k));
          f_dw.add_vertex(vertex_at(i + 1, j, k));
          f_dw.add_vertex(vertex_at(i + 1, j + 1, k));
        }

    for(int k = 0; k < nz; ++k)
      for(int j = 0; j < (periodic ? ny : ny + 1); ++j)
        for(int i = 0; i < nx; ++i) {
          //     0
          // 0 ----- 1
          //     1
          auto& e = edge_at(i, j, k, edge_color::horizontal);
          e = Edge(&e - edges_.data(), edge_color::horizontal);
          e.add_vertex(vertex_at(i, j, k));
          e.add_vertex(vertex_at(i + 1, j, k));

          if(j > 0 || periodic)
            e.add_face(face_at(i, j - 1, k, face_color::upward));
          if(j < ny || periodic)
            e.add_face(face_at(i, j, k, face_color::downward));
        }
    for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
        for(int i = 0; i < nx; ++i) {
          // 0
          //  \  0
          //   \ 
        // 1  \ 
        //     1
          auto& e = edge_at(i, j, k, edge_color::diagonal);
          e = Edge(&e - edges_.data(), edge_color::diagonal);
          e.add_vertex(vertex_at(i, j, k));
          e.add_vertex(vertex_at(i + 1, j + 1, k));

          e.add_face(face_at(i, j, k, face_color::downward));
          e.add_face(face_at(i, j, k, face_color::upward));
        }
    for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
        for(int i = 0; i < (periodic ? nx : nx + 1); ++i) {
          //     0
          // \^^^|\ 
        //  \1 | \ 
        //   \ | 0\ 
        //    \|___\ 
        //     1
          auto& e = edge_at(i, j, k, edge_color::vertical);
          e = Edge(&e - edges_.data(), edge_color::vertical);
          e.add_vertex(vertex_at(i, j, k));
          e.add_vertex(vertex_at(i, j + 1, k));

          if(i < nx || periodic)
            e.add_face(face_at(i, j, k, face_color::upward));
          if(i > 0 || periodic)
            e.add_face(face_at(i - 1, j, k, face_color::downward));
        }

    for(int k = 0; k < nz; ++k)
      for(int j = 0; j < (periodic ? ny : ny + 1); ++j)
        for(int i = 0; i < (periodic ? nx : nx + 1); ++i) {
          auto& v = vertex_at(i, j, k);
          //  1   2
          //   \  |
          //    \ |
          //     \|
          // 0---------- 3
          //      |\ 
        //      | \ 
        //      |  \ 
        //      5   4
          if(i > 0 || periodic) //
            v.add_edge(edge_at(i - 1, j, k, edge_color::horizontal));
          if(i > 0 && j > 0 || periodic) //
            v.add_edge(edge_at(i - 1, j - 1, k, edge_color::diagonal));
          if(j > 0 || periodic) //
            v.add_edge(edge_at(i, j - 1, k, edge_color::vertical));
          if(i < nx || periodic) //
            v.add_edge(edge_at(i, j, k, edge_color::horizontal));
          if(i < nx && j < ny || periodic) //
            v.add_edge(edge_at(i, j, k, edge_color::diagonal));
          if(j < ny || periodic) //
            v.add_edge(edge_at(i, j, k, edge_color::vertical));

          //    1
          //   \  |
          // 0  \ |  2
          //     \|
          //  ----------
          //      |\ 
          //   5  | \ 3
          //      |  \ 
          //        4
          if(i > 0 && j > 0 || periodic) {
            v.add_face(face_at(i - 1, j - 1, k, face_color::upward));
            v.add_face(face_at(i - 1, j - 1, k, face_color::downward));
          }
          if(i < nx && j > 0 || periodic) //
            v.add_face(face_at(i, j - 1, k, face_color::upward));
          if(i < nx && j < ny || periodic) {
            v.add_face(face_at(i, j, k, face_color::downward));
            v.add_face(face_at(i, j, k, face_color::upward));
          }
          if(i > 0 && j < ny || periodic)
            v.add_face(face_at(i - 1, j, k, face_color::downward));
        }
  }

  std::vector<Face> const& faces() const { return faces_; }
  std::vector<Vertex> const& vertices() const { return vertices_; }
  std::vector<Edge> const& edges() const { return edges_; }

  auto nx() const { return nx_; }
  auto ny() const { return ny_; }
  auto nz() const { return nz_; }

private:
  std::vector<Face> faces_;
  std::vector<Vertex> vertices_;
  std::vector<Edge> edges_;

  int nx_;
  int ny_;
  int nz_;
};

template <typename O, typename T>
class Data {
public:
  explicit Data(size_t size) : data_(size) {}
  T& operator()(O const& f) { return data_[f.id()]; }
  T const& operator()(O const& f) const { return data_[f.id()]; }
  auto begin() { return data_.begin(); }
  auto end() { return data_.end(); }

private:
  std::vector<T> data_;
};
template <typename T>
class FaceData : public Data<Face, T> {
public:
  explicit FaceData(Grid const& grid) : Data<Face, T>(grid.faces().size()) {}
};
template <typename T>
class VertexData : public Data<Vertex, T> {
public:
  explicit VertexData(Grid const& grid) : Data<Vertex, T>(grid.vertices().size()) {}
};
template <typename T>
class EdgeData : public Data<Edge, T> {
public:
  explicit EdgeData(Grid const& grid) : Data<Edge, T>(grid.edges().size()) {}
};
namespace faces {
template <class T, class Map, class Reduce>
auto reduce_on_vertices(VertexData<T> const& v_data, Grid const& grid, Map const& map,
                        Reduce const& reduce, FaceData<T> ret) {
  for(auto f : grid.faces()) {
    for(auto v : f.vertices())
      ret(f) = reduce(ret(f), map(v_data(*v)));
    ret(f) /= f.vertices().size();
  }
  return ret;
}
inline int edge_sign(Edge const& curr, Edge const& prev) {
  auto p0 = prev.vertex(0).id();
  auto p1 = prev.vertex(1).id();

  auto c0 = curr.vertex(0).id();
  auto c1 = curr.vertex(1).id();

  auto start = c0 == p0 || c0 == p1 ? c0 : c1;
  auto end = c0 == p0 || c0 == p1 ? c1 : c0;
  return start < end ? 1 : -1;
}
template <class T, class Map, class Reduce>
std::enable_if_t<(sizeof(std::declval<Map>()(std::declval<T>(), 1)) > 0), FaceData<float>>
reduce_on_edges(EdgeData<T> const& e_data, Grid const& grid, Map const& map, Reduce const& reduce,
                FaceData<float> ret) {
  for(auto f : grid.faces()) {
    // inward, if edge from negative to positive
    //   .
    //   |\ 
    //  0| \ 1   0: outwards
    //   |  \    1: inwards
    //   ----'   2: outwards
    //     2
    //
    //     1
    //   ----
    //   \  |
    //  0 \ |2   0: outwards
    //     \|    1: inwards
    //      ^    2: inwards
    auto edges = f.edges();
    for(auto prev = edges.end() - 1, curr = edges.begin(); curr != edges.end(); prev = curr, ++curr)
      ret(f) = reduce(ret(f), map(e_data(**curr), edge_sign(**curr, **prev)));
    ret(f) /= f.edges().size();
  }
  return ret;
}
template <class T, class Map, class Reduce>
std::enable_if_t<(sizeof(std::declval<Map>()(std::declval<T>())) > 0), FaceData<T>>
reduce_on_edges(EdgeData<T> const& e_data, Grid const& grid, Map const& map, Reduce const& reduce,
                FaceData<T> ret) {
  for(auto f : grid.faces()) {
    for(auto e : f.edges())
      ret(f) = reduce(ret(f), map(e_data(*e)));
    ret(f) /= f.edges().size();
  }
  return ret;
}

template <class T, class Map, class Reduce>
auto reduce_on_faces(FaceData<T> const& f_data, Grid const& grid, Map const& map,
                     Reduce const& reduce, FaceData<T> ret) {
  for(auto f : grid.faces()) {
    for(auto next_f : f.faces())
      ret(f) = reduce(ret(f), map(f_data(*next_f)));
    ret(f) /= f.faces().size();
  }
  return ret;
}

} // namespace faces

namespace edges {
template <class T, class Map, class Reduce>
auto reduce_on_faces(FaceData<T> const& f_data, Grid const& grid, Map const& map,
                     Reduce const& reduce, EdgeData<T> ret) {
  for(auto e : grid.edges()) {
    if(e.faces().size() == 2) {
      ret(e) = reduce(ret(e), map(f_data(e.face(0)), 1));
      ret(e) = reduce(ret(e), map(f_data(e.face(1)), -1));
    } else
      ret(e) = 0;
  }
  return ret;
}
} // namespace edges
std::ostream& toVtk(Grid const& grid, std::ostream& os = std::cout);
std::ostream& toVtk(std::string const& name, FaceData<double> const& f_data, Grid const& grid,
                    std::ostream& os = std::cout);
std::ostream& toVtk(std::string const& name, EdgeData<double> const& e_data, Grid const& grid,
                    std::ostream& os = std::cout);
std::ostream& toVtk(std::string const& name, VertexData<double> const& v_data, Grid const& grid,
                    std::ostream& os = std::cout);

namespace {
bool inner_face(Face const& f) {
  return (f.color() == face_color::downward && f.vertex(0).id() < f.vertex(1).id() &&
          f.vertex(0).id() < f.vertex(2).id()) ||
         (f.color() == face_color::upward && f.vertex(1).id() > f.vertex(0).id() &&
          f.vertex(1).id() > f.vertex(2).id());
}
} // namespace

} // namespace mylib
