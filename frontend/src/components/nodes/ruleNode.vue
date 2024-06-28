<template>
    <div>
        <div class="nodeBox" ref="nodeBox">
            <input type="text" v-model="nodeTitle" df-title>
        </div>
    </div>
</template>

<script>
export default {
    data() {
        return {
            nodeId: null,
            nodeTitle: '',
        };
    },
    mounted() {
    },
};
</script>

<style scoped>
.nodeBox {
    width: 100%;
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    text-align: center;
}

.nodeBox input {
    width: 100%;
    height: 100%;
    border: none;
    background-color: transparent;
    text-align: center;
    font-size: 1rem;
    font-weight: bold;
    color: #000;
}

.nodeBox input:focus {
    outline: none;
}
</style>