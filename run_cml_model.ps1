$workflow = if ($args.Count -ge 1) { $args[0] } else { "figures" }

docker build -t cml-model:paper .

docker run --rm `
  --env JULIA_NUM_THREADS=1 `
  -v "${PWD}/output:/app/output" `
  cml-model:paper $workflow
