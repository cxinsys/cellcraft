<template>
    <div>
        <div class="nodeBox" ref="nodeBox" @dblclick="startEditing">
            <input
                type="text"
                v-model="nodeTitle"
                ref="titleInput"
                @blur="stopEditing"
                @keyup.enter="stopEditing"
                @focus="startEditing"
                :readonly="!isEditing"
                df-title
            >
        </div>
    </div>
</template>

<script>
export default {
    data() {
        return {
            nodeId: null,
            nodeTitle: '',
            isEditing: false,
        };
    },
    methods: {
        startEditing() {
            this.isEditing = true;
            this.$nextTick(() => {
                if (this.$refs.titleInput) {
                    this.$refs.titleInput.focus();
                    this.$refs.titleInput.select();
                }
            });
        },
        stopEditing() {
            this.isEditing = false;
        },
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
    cursor: move;
    box-sizing: border-box;
}

.nodeBox input {
    width: 100%;
    height: 100%;
    border: none;
    background-color: transparent;
    text-align: center;
    font-size: 0.85rem;
    font-weight: bold;
    color: #333;
    padding: 2px 5px;
    box-sizing: border-box;
    cursor: move;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
}

.nodeBox input:readonly {
    cursor: move;
    user-select: none;
}

.nodeBox input:not(:readonly) {
    border: 1px solid #4a90e2;
    background-color: white;
    cursor: text;
    border-radius: 4px;
    font-size: 1rem;
    padding: 0 10px;
}

.nodeBox input:focus {
    outline: none;
    box-shadow: 0 0 5px #4a90e2;
}
</style>