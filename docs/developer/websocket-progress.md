# MetaPilot Web Progress Status

This note records the current web progress/API implementation for `1.1.0`.

## Current Role

`metapilot-web` is a Spring Boot wrapper that exposes a browser UI, REST endpoints and WebSocket/STOMP progress messages. It launches a MetaPilot JAR as a subprocess from `WorkflowService`.

For production Linux/HPC execution, prefer `metapilot-cli` directly. The web module is useful as an API/UI front end, but it is not the canonical scheduler runner.

## Implemented Components

| Component | Path | Status |
| --- | --- | --- |
| Spring Boot app | `metapilot-web/src/main/java/uk/ac/quadram/metapilot/web/MetaPilotWebApp.java` | Present |
| REST job controller | `web/rest/JobController.java` | Present |
| File upload controller | `web/rest/FileUploadController.java` | Present |
| Workflow service | `web/service/WorkflowService.java` | Runs a JAR subprocess |
| WebSocket config | `web/config/WebSocketConfig.java` | STOMP/SockJS endpoint |
| Progress broadcast service | `web/websocket/ProgressBroadcastService.java` | Sends progress/completion/cancellation messages |
| Static UI | `src/main/resources/static/index.html` | Present |
| Parameter builder UI | `src/main/resources/static/parameter-builder.html` | Present |

## Progress Message Shape

```json
{
  "jobId": "550e8400-e29b-41d4-a716-446655440000",
  "status": "RUNNING",
  "progress": 0.45,
  "message": "Processing",
  "timestamp": "2026-05-17T10:30:00"
}
```

Statuses include:

- `RUNNING`
- `COMPLETED`
- `FAILED`
- `CANCELLED`

## Endpoints

Typical job endpoints:

```text
GET  /api/jobs/status
POST /api/jobs/submit
GET  /api/jobs/{jobId}/status
POST /api/jobs/{jobId}/cancel
```

WebSocket/STOMP:

```text
/ws-stomp
/topic/jobs/{jobId}/progress
```

## Current Execution Caveat

`WorkflowService` locates a MetaPilot JAR and starts:

```text
java -jar <jar> --params <paramFile> --resources <resourceFile>
```

That wrapper path must be kept aligned with the CLI command surface. The current CLI supports subcommands such as `run`, `params-validate`, `dia-marker` and `mag-fragpipe`; the web wrapper should be updated before production use if it is expected to call the current `metapilot-cli` JAR.

## Production Hardening Still Needed

- Restrict WebSocket origins.
- Add authentication and job ownership checks.
- Persist job history if users need to reconnect after server restart.
- Route long-running work through a scheduler or queue rather than inside the web JVM.
- Update subprocess invocation to use the current CLI command contract.
- Avoid relying on CDN-hosted JavaScript assets for offline/HPC networks.
- Add integration tests for submit, progress, cancellation and failure paths.

## Build

```bash
mvn -pl metapilot-web -am -DskipTests package
```

Run locally:

```bash
java -jar metapilot-web/target/metapilot-web-1.1.0.jar
```

Then open:

```text
http://localhost:8080
```
