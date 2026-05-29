package uk.ac.quadram.metapilot.web.config;

import org.springframework.context.annotation.Configuration;
import org.springframework.messaging.simp.config.MessageBrokerRegistry;
import org.springframework.web.socket.config.annotation.EnableWebSocketMessageBroker;
import org.springframework.web.socket.config.annotation.StompEndpointRegistry;
import org.springframework.web.socket.config.annotation.WebSocketMessageBrokerConfigurer;

/**
 * WebSocket configuration for real-time progress updates.
 * 
 * <p>This configuration enables STOMP over WebSocket for bidirectional
 * communication between the server and web clients. Clients can subscribe
 * to progress updates for their jobs using the /topic/jobs/{jobId}/progress
 * destination.</p>
 * 
 * @author MetaPilot Team
 * @version 1.0.0
 * @since 1.0.0
 */
@Configuration
@EnableWebSocketMessageBroker
public class WebSocketConfig implements WebSocketMessageBrokerConfigurer {
    
    /**
     * Configure the message broker for WebSocket messaging.
     * 
     * <p>Destinations:
     * - /topic: For broadcasting messages to all subscribers
     * - /queue: For point-to-point messaging
     * </p>
     */
    @Override
    public void configureMessageBroker(MessageBrokerRegistry config) {
        // Enable a simple in-memory broker for destinations prefixed with "/topic"
        config.enableSimpleBroker("/topic", "/queue");
        
        // Prefix for messages bound to methods annotated with @MessageMapping
        config.setApplicationDestinationPrefixes("/app");
    }
    
    /**
     * Register STOMP endpoints.
     * 
     * <p>Clients connect to /ws-stomp endpoint with SockJS fallback.
     * This allows connections from any origin for development.</p>
     */
    @Override
    public void registerStompEndpoints(StompEndpointRegistry registry) {
        // Register /ws-stomp endpoint with SockJS fallback
        registry.addEndpoint("/ws-stomp")
                .setAllowedOriginPatterns("*") // Allow all origins for development
                .withSockJS();
    }
}

