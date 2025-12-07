#!/bin/bash

# CellCraft Installation Check Script
# Compact validation for CellCraft system setup

# Color codes
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly CYAN='\033[0;36m'
readonly NC='\033[0m'

# Configuration
readonly SERVICES=("frontend" "backend" "db" "rabbitmq" "celery")
readonly CPU_PLUGINS=("GENIE3" "GRNBoost2" "GRNViz" "LEAP" "Scribe" "TENET")
readonly GPU_PLUGINS=("GENIE3" "GRNBoost2" "GRNViz" "LEAP" "Scribe" "TENET" "FastSCODE" "FastTENET")
readonly TIMEOUT=10

# Runtime variables
MODE=""
PASSED=0
FAILED=0

# Output functions
pass() { echo -e "${GREEN}[PASS]${NC} $1"; ((PASSED++)); }
fail() { echo -e "${RED}[FAIL]${NC} $1"; ((FAILED++)); }
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }

# Print header
print_header() {
    echo -e "${CYAN}"
    echo "================================================"
    echo "        CellCraft Installation Check"
    echo "================================================"
    echo -e "${NC}"
}

# Mode selection
select_mode() {
    echo -e "${BLUE}Select installation mode:${NC}"
    echo "  1) CPU mode (6 plugins)"
    echo "  2) GPU mode (8 plugins)"
    echo ""
    read -p "Enter choice [1-2]: " choice

    case "$choice" in
        1) MODE="cpu" ;;
        2) MODE="gpu" ;;
        *)
            echo -e "${RED}Invalid choice. Exiting.${NC}"
            exit 1
            ;;
    esac

    info "Checking ${MODE^^} mode installation"
}

# Check container status
check_containers() {
    info "Checking container status..."

    local containers
    containers=$(docker ps --format '{{.Names}}' 2>/dev/null | grep cellcraft)

    if [[ -z "$containers" ]]; then
        fail "No CellCraft containers found running"
        return
    fi

    for service in "${SERVICES[@]}"; do
        if echo "$containers" | grep -q "$service"; then
            pass "${service} container running"
        else
            fail "${service} container not found"
        fi
    done
}

# Test connectivity
check_connectivity() {
    info "Testing HTTP endpoints..."

    # Frontend (8080)
    if curl -sf --connect-timeout "$TIMEOUT" http://localhost:8080 &>/dev/null; then
        pass "Frontend accessible (port 8080)"
    else
        fail "Frontend not accessible (port 8080)"
    fi

    # Backend API (8000)
    if curl -sf --connect-timeout "$TIMEOUT" http://localhost:8000/docs &>/dev/null; then
        pass "Backend API accessible (port 8000)"
    else
        fail "Backend API not accessible (port 8000)"
    fi
}

# Verify plugins
verify_plugins() {
    info "Verifying plugins for ${MODE^^} mode..."

    local expected_plugins
    if [[ "$MODE" == "cpu" ]]; then
        expected_plugins=("${CPU_PLUGINS[@]}")
    else
        expected_plugins=("${GPU_PLUGINS[@]}")
    fi

    local response
    response=$(curl -sf --connect-timeout "$TIMEOUT" http://localhost:8000/api/plugin/list 2>/dev/null)

    if [[ -z "$response" ]]; then
        fail "Cannot access plugin registry"
        return
    fi

    local found_count=0
    local missing_plugins=()

    for plugin in "${expected_plugins[@]}"; do
        if echo "$response" | grep -qi "\"name\"[[:space:]]*:[[:space:]]*\"${plugin}\""; then
            ((found_count++))
        else
            missing_plugins+=("$plugin")
        fi
    done

    if [[ $found_count -eq ${#expected_plugins[@]} ]]; then
        pass "All ${#expected_plugins[@]} plugins registered"
    else
        fail "Found $found_count/${#expected_plugins[@]} expected plugins"
        if [[ ${#missing_plugins[@]} -gt 0 ]]; then
            warn "Missing plugins: ${missing_plugins[*]}"
        fi
    fi
}

# Print summary
print_summary() {
    echo ""
    echo -e "${CYAN}================================================${NC}"
    echo -e "${CYAN}                    SUMMARY${NC}"
    echo -e "${CYAN}================================================${NC}"
    echo -e "Mode:   ${CYAN}${MODE^^}${NC}"
    echo -e "Passed: ${GREEN}${PASSED}${NC}"
    echo -e "Failed: ${RED}${FAILED}${NC}"
    echo ""

    if [[ $FAILED -eq 0 ]]; then
        echo -e "${GREEN}Installation verified successfully!${NC}"
        echo ""
        echo -e "${BLUE}Access Points:${NC}"
        echo "  Frontend:    http://localhost:8080"
        echo "  Backend API: http://localhost:8000/docs"
    else
        echo -e "${RED}Installation has issues. Review failed checks above.${NC}"
        echo ""
        echo -e "${BLUE}Troubleshooting:${NC}"
        echo "  Docs:    https://cellcraft.gitbook.io/cellcraft-docs"
        echo "  Logs:    docker compose logs"
        echo "  Restart: docker compose restart"
    fi
    echo ""
}

# Main function
main() {
    print_header
    select_mode
    echo ""
    check_containers
    echo ""
    check_connectivity
    echo ""
    verify_plugins
    print_summary

    if [[ $FAILED -eq 0 ]]; then
        exit 0
    else
        exit 1
    fi
}

main "$@"
