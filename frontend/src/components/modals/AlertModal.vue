<template>
    <div v-if="show" class="modal-overlay">
        <div class="modal-content">
            <p class="modal-title">{{ title }}</p>
            <ul>
                <li v-for="(message, index) in messages" :key="index">{{ message }}</li>
            </ul>
            <button @click="closeAlertModal">Confirm</button>
        </div>
    </div>
</template>

<script>
export default {
    props: {
        show: {
            type: Boolean,
            required: true
        },
        title: {
            type: String,
            required: true
        },
        messages: {
            type: Array,
            required: true
        }
    },
    methods: {
        closeAlertModal() {
            this.$emit('close');
        }
    }
};
</script>

<style scoped>
.modal-overlay {
    position: fixed;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    background-color: rgba(0, 0, 0, 0.5);
    display: flex;
    justify-content: center;
    align-items: center;
    z-index: 9999;
}

.modal-content {
    background-color: #fff;
    padding: 1rem;
    border-radius: 0.8rem;
    text-align: left;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
}

.modal-title {
    font-size: 1.1rem;
    font-weight: bold;
    margin-bottom: 1rem;
    color: red;
}

ul {
    list-style-type: disc;
    margin-left: 1rem;
    margin-bottom: 1rem;
    display: flex;
    flex-direction: column;
    justify-content: center;
    gap: 0.5rem;
}

ul li {
    font-size: 1rem;
    line-height: 1.1rem;
}

button {
    background-color: #007BFF;
    color: white;
    padding: 0.5rem 1rem;
    border: none;
    border-radius: 0.5rem;
    cursor: pointer;
}

button:hover {
    background-color: #0056b3;
}
</style>