function loops = find_independent_loops(measurements)
  % Finds independent loops in a network of levelling measurements.

  % Args:
  %   measurements: A list of levelling measurements, each represented as a
  %     tuple of (A, B, dH, dist, n, diff).

  % Returns:
  %   A list of independent loops. Each loop is a list of levelling measurements.

  loops = [];
  visited = set();
  for measurement in measurements:
    if measurement[0] not in visited:
      loop = [measurement];
      visited.add(measurement[0]);
      for measurement in measurements:
        if measurement[1] == loop[-1][0] and measurement[0] not in visited:
          loop.append(measurement);
          visited.add(measurement[0]);
      loops.append(loop);
  return loops;
end

if __name__ == "__main__":
  measurements = [
      ("A", "B", 10, 100, 2, 0),
      ("B", "C", 5, 50, 2, 0),
      ("C", "D", 15, 150, 2, 0),
      ("D", "A", 10, 100, 2, 0),
      ("A", "E", 20, 200, 2, 0),
      ("E", "F", 10, 100, 2, 0),
      ("F", "G", 5, 50, 2, 0),
      ("G", "A", 15, 150, 2, 0),
  ]
  loops = find_independent_loops(measurements);
  print(loops);

end
