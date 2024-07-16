<template>
    <div class="control-popup__files" v-if="show_files">
        <table class="control-popup__table">
            <thead>
                <tr>
                    <th>Name</th>
                    <th>Date</th>
                    <th>Type</th>
                    <th>Size</th>
                </tr>
            </thead>
            <tbody>
                <tr v-for="(file, index) in files_list" :key="index">
                    <td>{{ file.file_name | cutFromDotName }}</td>
                    <td>{{ file.created_at | cutFromT }}</td>
                    <td>{{ file.file_name | cutFromDotType }}</td>
                    <td>{{ file.file_size | formatBytes }}</td>
                </tr>
            </tbody>
        </table>
    </div>
</template>

<script>
export default {
    props: {
        show_files: {
            type: Boolean,
            required: true
        },
        files_list: {
            type: Array,
            required: true
        }
    },
    filters: {
        cutFromDotName(value) {
            // Implement the filter logic here
            return value.split('.')[0];
        },
        cutFromT(value) {
            // Implement the filter logic here
            return value.split('T')[0];
        },
        cutFromDotType(value) {
            // Implement the filter logic here
            return value.split('.').pop();
        },
        formatBytes(value) {
            // Implement the filter logic here
            if (value === 0) return '0 Bytes';
            const k = 1024;
            const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
            const i = Math.floor(Math.log(value) / Math.log(k));
            return parseFloat((value / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
        }
    }
};
</script>

<style scoped>
.control-popup__files {
    /* width: 8rem; */
    width: 40vw;
    max-width: 400px;
    /* height: 34rem; */
    height: 30vh;
    max-height: 300px;

    border-radius: 16px;
    background: rgba(244, 246, 251, 0.586);
    box-shadow: 0px 0px 5px 0px rgba(0, 0, 0, 1);
    position: absolute;
    bottom: 98px;
    right: calc(50% + 1vw);
    z-index: 9998;
    opacity: 1;
    display: flex;
    align-items: center;
    justify-content: center;
}

.control-popup__table {
  width: 95%;
  height: auto;
  margin: auto;
  border-collapse: collapse;
  position: absolute;
  top: 20px;
}

.control-popup__table thead {
  height: 26px;
  font-weight: 500;
  color: rgb(49, 49, 49);
  border-bottom: 1px solid #6767678c;
}

.control-popup__table td {
  vertical-align: middle;
  font-weight: 400;
  text-align: center;
  color: rgb(68, 68, 68);
  padding: 0.7rem;
  margin: 1rem;
}

.control-popup__files::-webkit-scrollbar {
  width: 10px;
  /* width of the entire scrollbar */
}

.control-popup__files::-webkit-scrollbar-track {
  background: #f1f1f1;
  /* color of the tracking area */
  border-radius: 16px;
  /* keep the same radius as the container */
}

.control-popup__files::-webkit-scrollbar-thumb {
  background: #888;
  /* color of the scroll thumb */
  border-radius: 16px;
  /* keep the same radius as the container */
}

.control-popup__files::-webkit-scrollbar-thumb:hover {
  background: #555;
  /* color of the scroll thumb on hover */
}
</style>