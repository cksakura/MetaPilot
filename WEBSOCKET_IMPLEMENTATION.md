# WebSocket Real-Time Progress Implementation

## Overview

MetaPilot Web now includes WebSocket support for real-time progress updates during workflow execution. Users can see live progress bars and status messages without manually refreshing the page.

## Architecture

### Server-Side (Spring Boot)

1. **WebSocketConfig** (`config/WebSocketConfig.java`)
   - Configures STOMP messaging over WebSocket
   - Enables message broker for `/topic` and `/queue` destinations
   - Registers `/ws-stomp` endpoint with SockJS fallback

2. **ProgressBroadcastService** (`websocket/ProgressBroadcastService.java`)
   - Service for sending progress updates via WebSocket
   - Uses Spring's `SimpMessagingTemplate` to broadcast messages
   - Supports progress updates, completions, and cancellations

3. **ProgressMessage** (`websocket/ProgressMessage.java`)
   - Data class for progress updates
   - Contains: jobId, status, progress (0.0-1.0), message, timestamp

4. **WorkflowService Integration**
   - Injected `ProgressBroadcastService` into `WorkflowService`
   - Broadcasts progress at key workflow stages
   - Handles cancellation and completion updates

### Client-Side (JavaScript)

1. **WebSocket Connection**
   - Uses SockJS + STOMP library
   - Connects to `/ws-stomp` endpoint
   - Auto-reconnects on connection loss

2. **Progress Subscriptions**
   - Subscribes to `/topic/jobs/{jobId}/progress` for each job
   - Receives `ProgressMessage` objects
   - Updates UI in real-time

3. **UI Updates**
   - Progress bars update smoothly
   - Status badges change color based on state
   - Messages display current workflow stage
   - Auto-removal of completed jobs after 10 seconds

## Message Types

### Progress Update
```json
{
  "jobId": "550e8400-e29b-41d4-a716-446655440000",
  "status": "RUNNING",
  "progress": 0.45,
  "message": "Processing... 45%",
  "timestamp": "2024-01-15T10:30:00"
}
```

### Completion
```json
{
  "jobId": "550e8400-e29b-41d4-a716-446655440000",
  "status": "COMPLETED",
  "progress": 1.0,
  "message": "Job completed successfully",
  "timestamp": "2024-01-15T10:35:00"
}
```

### Failure
```json
{
  "jobId": "550e8400-e29b-41d4-a716-446655440000",
  "status": "FAILED",
  "progress": 0.0,
  "message": "Job failed: Error message",
  "timestamp": "2024-01-15T10:32:00"
}
```

### Cancellation
```json
{
  "jobId": "550e8400-e29b-41d4-a716-446655440000",
  "status": "CANCELLED",
  "progress": 0.0,
  "message": "Job cancelled by user",
  "timestamp": "2024-01-15T10:31:00"
}
```

## Integration Points

### WorkflowService

The `WorkflowService` broadcasts progress at these stages:

1. **Job Start**: `broadcastProgress(jobId, "RUNNING", 0.0, "Job started")`
2. **During Execution**: Periodic updates with progress percentage
3. **Job Completion**: `broadcastCompletion(jobId, true, "Job completed successfully")`
4. **Job Failure**: `broadcastCompletion(jobId, false, "Job failed: error")`
5. **Job Cancellation**: `broadcastCancellation(jobId, "Job cancelled by user")`

### Frontend Connection

```javascript
// Connect to WebSocket
function connectWebSocket() {
    const socket = new SockJS(window.location.origin + '/ws-stomp');
    stompClient = Stomp.over(socket);
    
    stompClient.connect({}, function(frame) {
        console.log('WebSocket connected: ' + frame);
    });
}

// Subscribe to job progress
function subscribeToJob(jobId) {
    stompClient.subscribe('/topic/jobs/' + jobId + '/progress', function(message) {
        const progress = JSON.parse(message.body);
        updateJobProgress(jobId, progress);
    });
}
```

## Features

### ✅ Implemented
- Real-time progress bars
- Live status updates
- Animated progress transitions
- Status badges (running/completed/failed/cancelled)
- Auto-removal of completed jobs
- WebSocket reconnection handling
- Multiple concurrent job tracking

### 🔄 UI Behavior
- Progress bars smoothly animate from 0% to 100%
- Status color coding:
  - **Running**: Yellow/Orange
  - **Completed**: Green
  - **Failed/Cancelled**: Red
- Jobs auto-remove 10 seconds after completion
- "No jobs" message appears when all jobs complete
- Cancel buttons only visible for running jobs

## Testing

### Manual Testing

1. Start the server:
```bash
java -jar metapilot-web/target/metapilot-web-1.0.0-beta.2.jar
```

2. Open browser console (F12)

3. Navigate to: `http://localhost:8080`

4. Submit a test job:
   - Parameter File: `/test/params.json`
   - Resource File: `/test/resources.params`
   - Click "Submit Job"

5. Observe:
   - WebSocket connection in console
   - Real-time progress updates
   - Progress bar animation
   - Status changes

6. Cancel a job:
   - Click "Cancel" button
   - Observe cancellation message
   - Job status changes to CANCELLED

### Expected Console Output

```
WebSocket connected: CONNECTED
```

## Configuration

### Server Configuration

Edit `application.properties`:

```properties
# WebSocket is enabled by default
# No additional configuration needed
```

### Client Configuration

The frontend uses CDN-hosted libraries:
- SockJS 1.6.1
- STOMP 2.3.3

To use local libraries, download and serve from `/static/js/`.

## Troubleshooting

### WebSocket Not Connecting

**Symptom**: No console log for "WebSocket connected"

**Possible Causes**:
1. Browser blocked WebSocket (check firewall)
2. SockJS CDN blocked (check network tab)
3. CORS issue (already configured to allow all origins)

**Solution**: Check browser console for errors, verify SockJS CDN is accessible

### Progress Not Updating

**Symptom**: Progress bar stuck at 0%

**Possible Causes**:
1. Not subscribed to job topic
2. Broadcast service not injected
3. WebSocket disconnected

**Solution**: Check WebSocket connection status in console, verify job ID matches subscription

### Multiple Concurrent Jobs

**Symptom**: Only one job shows progress

**Solution**: Each job requires separate subscription. The code automatically subscribes when jobs are submitted.

## Future Enhancements

- [ ] Job history persistence
- [ ] Download results link
- [ ] Estimated time remaining
- [ ] Detailed log viewer
- [ ] Email notifications on completion
- [ ] WebSocket authentication
- [ ] Compression for large messages

## Dependencies

### Backend
- Spring Boot WebSocket (`spring-boot-starter-websocket`)
- SockJS support built-in

### Frontend
- SockJS Client 1.6.1 (CDN)
- STOMP.js 2.3.3 (CDN)

## Performance

- Message broker: In-memory (sufficient for single-instance deployment)
- Broadcast frequency: ~1 second per update (configurable)
- Memory usage: Minimal (messages are small JSON)
- Thread safety: All operations are thread-safe

## Security Considerations

⚠️ **Current Setup**: Allows connections from any origin

**For Production**:
1. Restrict allowed origins in `WebSocketConfig`
2. Add authentication/authorization
3. Validate job ownership before broadcasting
4. Rate limit connections

## Example Usage Flow

1. User submits job → JobController creates job ID
2. WorkflowService starts execution in background thread
3. WorkflowService broadcasts "Job started" (0%)
4. During execution, broadcasts progress updates (10%, 20%, ...)
5. Frontend receives updates → Updates progress bar
6. On completion, broadcasts "Job completed" (100%)
7. Frontend hides job after 10 seconds

This creates a seamless real-time monitoring experience for users!

