# Benchmark Report

## Best configurations
- Strict (1e-5): depth=4, iterations=1500, seed=10, end_len=228, runtime_sec=29.8365\n
- Loose (1e-3): depth=4, iterations=1500, seed=5, end_len=226, runtime_sec=2.4003\n
## Combined table

| depth | iterations | runs | best_strict | mean_strict | eq_strict | best_loose | mean_loose | eq_loose |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 500 | 3 | 266 | 274 | 1.0 | 266 | 274 | 1.0 |
| 3 | 1000 | 3 | 246 | 250 | 0.667 | 246 | 250 | 1.0 |
| 3 | 1500 | 3 | 240 | 242.67 | 0.667 | 240 | 242.67 | 1.0 |
| 4 | 500 | 3 | 262 | 265.33 | 1.0 | 262 | 265.33 | 1.0 |
| 4 | 1000 | 3 | 238 | 243.33 | 1.0 | 238 | 243.33 | 1.0 |
| 4 | 1500 | 3 | 226 | 231.33 | 0.667 | 226 | 231.33 | 1.0 |
