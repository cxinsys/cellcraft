<template>
    <section class="control-bar">
        <ul class="control-bar__btnList">
            <li class="control-bar__button" @click="toggleFile">
                <img class="control-bar__icon" src="@/assets/control_files.png" />
            </li>
            <li class="control-bar__button" @click="downloadDrawflow">
                <img class="control-bar__icon" src="@/assets/control_save.png" />
            </li>
            <li class="control-bar__button">
                <button class="run_button" @click="runWorkflow">
                    <img class="control-bar__icon" src="@/assets/control_run.png" />
                </button>
            </li>
            <li>
                <div class="loader" @click="toggleTask" v-if="on_progress == true"></div>
                <div class="loader_done" @click="toggleTask" v-else></div>
            </li>
            <li class="control-bar__button">
                <img class="control-bar__icon white margin__top-4" v-if="isTabView" src="@/assets/view.png"
                    @click="toggleTabView" />
                <img class="control-bar__icon white margin__top-4" v-else src="@/assets/view_hide.png"
                    @click="toggleTabView" />
            </li>
        </ul>
    </section>
</template>

<script>
export default {
    props: {
        on_progress: {
            type: Boolean,
            required: true
        },
        isTabView: {
            type: Boolean,
            required: true
        }
    },
    methods: {
        toggleFile() {
            this.$emit('toggle-file');
        },
        saveWorkflowProject() {
            this.$emit('save-workflow-project');
        },
        runWorkflow() {
            this.$emit('run-workflow');
        },
        toggleTask() {
            this.$emit('toggle-task');
        },
        toggleTabView() {
            this.$emit('toggle-tab-view');
        },
        downloadDrawflow() {
            this.$emit('download-drawflow');
        }
    }
};
</script>

<style scoped>
.control-bar {
    height: 50px;
    width: 260px;
    border-radius: 10px;
    background: rgba(255, 255, 255, 0.1);
    box-shadow: 0px 0px 1px 0px rgba(255, 255, 255, 0.5);
    position: absolute;
    bottom: 24px;
    left: calc(50% - 130px);
    display: flex;
    align-items: center;
    justify-content: center;
}

.control-bar__btnList {
    width: 100%;
    height: 100%;
    display: flex;
    justify-content: center;
    align-items: center;
}

.control-bar__button {
    width: 24px;
    height: 24px;
    margin: 0 8px;
    align-items: center;
}

.control-bar__icon {
    max-width: 24px;
    max-height: 24px;
    object-fit: contain;
    object-position: center;
    opacity: 0.6;
}

.run_button {
  width: 100%;
  height: 100%;
  background: none;
  border: none;
  cursor: pointer;
}

.loader,
.loader_done {
  border: 4px solid #f3f3f3bf;
  border-radius: 50%;
  margin-left: 8px;
  margin-right: 6px;
  width: 20px;
  height: 20px;
  opacity: 0.5;
}

.loader {
  border-top: 4px solid #41b3ff;
  animation: spin 3s linear infinite;
}

.loader_done:hover {
  opacity: 0.7;
  /* You can adjust this value to your liking */
  transform: scale(1.1);
}
</style>