import heapq
import random

def dijkstra(graph, start, end):
    # Priority queue to store (distance, vertex) tuples
    queue = [(0, start)]
    # Dictionary to store the shortest path to each vertex
    distances = {vertex: float('infinity') for vertex in graph}
    distances[start] = 0
    # Dictionary to store the path
    previous_vertices = {vertex: None for vertex in graph}

    while queue:
        current_distance, current_vertex = heapq.heappop(queue)

        if current_distance > distances[current_vertex]:
            continue

        for neighbor, weight in graph[current_vertex].items():
            distance = current_distance + weight

            if distance < distances[neighbor]:
                distances[neighbor] = distance
                previous_vertices[neighbor] = current_vertex
                heapq.heappush(queue, (distance, neighbor))

    path, current_vertex = [], end
    while previous_vertices[current_vertex] is not None:
        path.insert(0, current_vertex)
        current_vertex = previous_vertices[current_vertex]
    if path:
        path.insert(0, current_vertex)
    return path

# Example usage:
# graph = {
#     'A': {'B': 1, 'C': 4},
#     'B': {'A': 1, 'C': 2, 'D': 5},
#     'C': {'A': 4, 'B': 2, 'D': 1},
#     'D': {'B': 5, 'C': 1}
# }

# start_vertex = 'A'
# end_vertex = 'D'
# print(dijkstra(graph, start_vertex, end_vertex))

def generate_random_graph(num_nodes, max_edges, min_cost=2, max_cost=20):
    graph = {i: {} for i in range(num_nodes)}
    edges = set()

    while len(edges) < max_edges:
        u = random.randint(0, num_nodes - 1)
        v = random.randint(0, num_nodes - 1)
        if u != v and (u, v) not in edges and (v, u) not in edges:
            cost = random.randint(min_cost, max_cost)
            graph[u][v] = cost
            graph[v][u] = cost
            edges.add((u, v))

    

    return graph, edges

def is_connected(graph):
    visited = set()
    def dfs(v):
        visited.add(v)
        for neighbor in graph[v]:
            if neighbor not in visited:
                dfs(neighbor)
    
    # Inicia a DFS a partir do primeiro vértice
    start_vertex = next(iter(graph))
    dfs(start_vertex)
    
    # Verifica se todos os vértices foram visitados
    return len(visited) == len(graph)


def find_path_and_write_to_file(num_nodes, max_edges, filename):
    while True:
        start_vertex = random.randint(0, num_nodes - 1)
        end_vertex = random.randint(0, num_nodes - 1)
        while start_vertex == end_vertex:
            end_vertex = random.randint(0, num_nodes - 1)

        graph,edges_set = generate_random_graph(num_nodes, max_edges)

        if not is_connected(graph):
            continue

        path = dijkstra(graph, start_vertex, end_vertex)
        print(graph)
        print(path)

        if path:
            with open(filename, 'w') as file:
                file.write(f"{num_nodes}\n")
                file.write(f"{start_vertex} {end_vertex}\n")
                file.write("\nEdges\n")

                for u in graph:
                    for v, cost in graph[u].items():
                        if u < v:  # To avoid writing both (u, v) and (v, u)
                            file.write(f"{u} {v} {cost}\n")

            # Cria um conjunto de arestas do caminho
            path_edges = set((path[i], path[i+1]) if path[i] < path[i+1] else (path[i+1], path[i]) for i in range(len(path) - 1))

            # Remove as arestas do caminho de edges_set
            edges_set -= path_edges

            toll_edges = random.sample(list(edges_set), min(max_toll_edges, len(edges_set)))
            with open(filename, 'a') as file:
                file.write("\nTollEdges\n")
                for u, v in toll_edges:
                    file.write(f"{u} {v}\n")
            break
        

number_of_nodes = 60
max_number_of_edges = 208
max_toll_edges = 42
#graph, edges_set = generate_random_graph(number_of_nodes,  max_number_of_edges)
# print(graph)
# print(edges_set)


find_path_and_write_to_file(number_of_nodes, max_number_of_edges, 'instancia20.txt')